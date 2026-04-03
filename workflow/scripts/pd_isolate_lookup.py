import csv
import io
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

import pandas as pd


FTP_RESULTS_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
FTP_SCANABLE_DISTANCE_BYTES = 2 * 1024 * 1024 * 1024
WORKING_COLUMNS = [
    "Sample",
    "srr_acc",
    "biosample_acc",
    "pd_target_acc",
    "pd_asm_acc",
    "scientific_name",
    "pd_snp_cluster",
    "pd_taxgroup",
    "pd_scientific_name",
    "pd_source_type",
    "pd_host",
    "pd_isolation_source",
    "pd_comparator_mode",
    "pd_comparator_count",
    "pd_lookup_status",
    "pd_lookup_source",
    "pd_lookup_note",
]

OUTPUT_COLUMNS = [
    "Sample",
    "srr_acc",
    "biosample_acc",
    "pd_target_acc",
    "pd_asm_acc",
    "pd_snp_cluster",
    "pd_taxgroup",
    "pd_scientific_name",
    "pd_source_type",
    "pd_host",
    "pd_isolation_source",
    "pd_comparator_mode",
    "pd_comparator_count",
    "pd_lookup_status",
    "pd_lookup_source",
    "pd_lookup_note",
]

COMPARATOR_COLUMNS = [
    "Sample",
    "pd_query_target_acc",
    "pd_snp_cluster",
    "pd_taxgroup",
    "comparator_rank",
    "comparator_mode",
    "comparator_distance",
    "comparator_target_acc",
    "comparator_biosample_acc",
    "comparator_asm_acc",
    "comparator_run",
    "comparator_sample_name",
    "comparator_scientific_name",
    "comparator_source_type",
    "comparator_host",
    "comparator_isolation_source",
]

SAMPLE_METADATA_ALIASES = {
    "Sample": ["sample", "sample_id", "sample name", "name"],
    "srr_acc": ["srr", "srr_acc", "run", "run_acc", "run accession", "sra", "sra_run"],
    "biosample_acc": ["biosample", "biosample_acc", "biosample accession", "biosample_accession"],
    "pd_target_acc": ["target_acc", "target accession", "target"],
    "pd_asm_acc": ["asm_acc", "assembly", "assembly accession", "assembly_acc"],
    "scientific_name": ["scientific_name", "scientific name", "organism", "species"],
}

PD_ISOLATE_ALIASES = {
    "biosample_acc": ["biosample_acc", "biosample", "biosample accession", "biosample_accession"],
    "pd_target_acc": ["target_acc", "target accession", "target"],
    "pd_asm_acc": ["asm_acc", "assembly", "assembly accession", "assembly_acc"],
    "pd_snp_cluster": ["erd_group", "snp cluster", "snp_cluster", "snp cluster id", "pds_acc"],
    "pd_taxgroup": ["taxgroup_name", "taxgroup", "organism", "taxgroup name"],
    "pd_scientific_name": ["scientific_name", "scientific name"],
    "pd_source_type": ["source_type", "source type"],
    "pd_host": ["host"],
    "pd_isolation_source": ["isolation_source", "isolation source"],
}

PD_EXCEPTION_ALIASES = {
    "biosample_acc": ["biosample_acc", "biosample", "biosample accession", "biosample_accession"],
    "srr_acc": ["run", "runs", "srr_acc", "srr", "run accession"],
    "pd_lookup_note": [
        "exception",
        "exception_reason",
        "failure_reason",
        "message",
        "reason",
        "description",
        "comment",
        "consequence",
    ],
}

ORGANISM_GROUP_OVERRIDES = {
    "escherichia coli": "Escherichia_coli_Shigella",
    "shigella boydii": "Escherichia_coli_Shigella",
    "shigella dysenteriae": "Escherichia_coli_Shigella",
    "shigella flexneri": "Escherichia_coli_Shigella",
    "shigella sonnei": "Escherichia_coli_Shigella",
}


def as_optional_path(value: object) -> Optional[Path]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"none", "null", "na"}:
        return None
    return Path(text)


def normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def match_columns(columns: list[str], aliases: dict[str, list[str]]) -> dict[str, str]:
    normalized = {normalize_header(col): col for col in columns}
    matched = {}
    for canonical, candidates in aliases.items():
        for candidate in candidates:
            original = normalized.get(normalize_header(candidate))
            if original is not None:
                matched[canonical] = original
                break
    return matched


def read_delimited_table(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        delimiter = dialect.delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample else ","
    return pd.read_csv(path, sep=delimiter, dtype=str, encoding="utf-8-sig").fillna("")


def infer_srr(value: str) -> str:
    return value if re.fullmatch(r"[SED]RR\d+", value) else ""


def infer_biosample(value: str) -> str:
    return value if re.fullmatch(r"SAM[A-Z0-9]+", value) else ""


def load_mash_species_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    mash = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "user_genome" not in mash.columns or "classification" not in mash.columns:
        return {}
    mash["Sample"] = mash["user_genome"].str.replace(".chromosome", "", regex=False)
    mash["scientific_name"] = mash["classification"].str.replace(r".*s__", "", regex=True)
    mash["scientific_name"] = mash["scientific_name"].replace("unknown", "")
    mash = mash[mash["Sample"] != ""]
    return dict(zip(mash["Sample"], mash["scientific_name"]))


def build_sample_records(samples: list[str], metadata_path: Optional[Path], mash_path: Path) -> pd.DataFrame:
    mash_species_map = load_mash_species_map(mash_path)
    records = pd.DataFrame({"Sample": samples}, dtype=str).assign(
        srr_acc=[infer_srr(sample) for sample in samples],
        biosample_acc=[infer_biosample(sample) for sample in samples],
        pd_target_acc="",
        pd_asm_acc="",
        scientific_name=[mash_species_map.get(sample, "") for sample in samples],
    )

    if metadata_path is None or not metadata_path.exists():
        return records

    metadata = read_delimited_table(metadata_path)
    column_map = match_columns(metadata.columns.tolist(), SAMPLE_METADATA_ALIASES)
    if "Sample" not in column_map:
        return records

    selected = pd.DataFrame({"Sample": metadata[column_map["Sample"]].astype(str)})
    for canonical in ["srr_acc", "biosample_acc", "pd_target_acc", "pd_asm_acc", "scientific_name"]:
        selected[canonical] = metadata[column_map[canonical]].astype(str) if canonical in column_map else ""

    merged = records.merge(selected, on="Sample", how="left", suffixes=("", "_meta")).fillna("")
    for canonical in ["srr_acc", "biosample_acc", "pd_target_acc", "pd_asm_acc", "scientific_name"]:
        meta_column = f"{canonical}_meta"
        merged[canonical] = merged[meta_column].where(merged[meta_column] != "", merged[canonical])

    merged["srr_acc"] = merged.apply(
        lambda row: row["srr_acc"] if row["srr_acc"] else infer_srr(row["Sample"]), axis=1
    )
    merged["biosample_acc"] = merged.apply(
        lambda row: row["biosample_acc"] if row["biosample_acc"] else infer_biosample(row["Sample"]),
        axis=1,
    )
    return merged[["Sample", "srr_acc", "biosample_acc", "pd_target_acc", "pd_asm_acc", "scientific_name"]]


def fetch_sra_metadata(accession: str) -> dict[str, str]:
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc={accession}"
    result = {"biosample_acc": "", "scientific_name": "", "taxid": "", "note": ""}
    try:
        with urlopen(url, timeout=20) as response:
            payload = response.read()
    except (HTTPError, URLError, TimeoutError) as exc:
        result["note"] = f"SRR lookup failed for {accession}: {exc}"
        return result

    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        result["note"] = f"SRR lookup returned malformed XML for {accession}: {exc}"
        return result

    biosample = root.findtext(".//EXTERNAL_ID[@namespace='BioSample']")
    scientific_name = root.findtext(".//SCIENTIFIC_NAME")
    taxid = root.findtext(".//TAXON_ID")
    result["biosample_acc"] = biosample.strip() if biosample else ""
    result["scientific_name"] = scientific_name.strip() if scientific_name else ""
    result["taxid"] = taxid.strip() if taxid else ""
    if not result["biosample_acc"]:
        result["note"] = f"No BioSample accession found in SRA metadata for {accession}"
    return result


def open_url(url: str, method: str = "GET"):
    request = Request(url, headers={"User-Agent": "SWAM-g PD lookup"}, method=method)
    return urlopen(request, timeout=60)


def fetch_directory_entries(url: str) -> list[str]:
    with open_url(url) as response:
        html = response.read().decode("utf-8", errors="replace")
    entries = re.findall(r'href="([^"]+)"', html)
    cleaned = []
    for entry in entries:
        name = entry.rstrip("/")
        if name in {"", ".."} or name.startswith("?") or name.startswith("/"):
            continue
        if name.startswith("https://www.hhs.gov/"):
            continue
        cleaned.append(name)
    return cleaned


def get_remote_size_bytes(url: str) -> Optional[int]:
    try:
        with open_url(url, method="HEAD") as response:
            value = response.headers.get("Content-Length")
    except (HTTPError, URLError, TimeoutError, ValueError):
        return None
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def stream_tsv_rows(url: str):
    with open_url(url) as response:
        with io.TextIOWrapper(response, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                cleaned = {}
                for key, value in row.items():
                    text = (value or "").strip()
                    cleaned[key] = "" if text.upper() == "NULL" else text
                yield cleaned


def load_pd_table(path: Optional[Path], aliases: dict[str, list[str]]) -> tuple[Optional[pd.DataFrame], str]:
    if path is None:
        return None, ""
    if not path.exists():
        return None, f"File not found: {path}"

    table = read_delimited_table(path)
    column_map = match_columns(table.columns.tolist(), aliases)
    if "biosample_acc" not in column_map and "srr_acc" not in column_map:
        return None, f"File missing accession column: {path}"

    renamed = pd.DataFrame(
        {
            "biosample_acc": table[column_map["biosample_acc"]].astype(str) if "biosample_acc" in column_map else "",
            "srr_acc": table[column_map["srr_acc"]].astype(str) if "srr_acc" in column_map else "",
        }
    )
    for canonical in aliases:
        if canonical in {"biosample_acc", "srr_acc"}:
            continue
        renamed[canonical] = table[column_map[canonical]].astype(str) if canonical in column_map else ""

    renamed["biosample_acc"] = renamed["biosample_acc"].str.strip()
    renamed["srr_acc"] = renamed["srr_acc"].str.strip()
    renamed = renamed[(renamed["biosample_acc"] != "") | (renamed["srr_acc"] != "")]
    renamed = renamed.groupby(["biosample_acc", "srr_acc"], as_index=False).first()
    return renamed.fillna(""), ""


class PDLiveClient:
    def __init__(self):
        self.organism_groups: Optional[list[str]] = None
        self.group_releases: dict[str, dict[str, str]] = {}

    def list_organism_groups(self) -> list[str]:
        if self.organism_groups is None:
            self.organism_groups = [entry for entry in fetch_directory_entries(FTP_RESULTS_URL) if "." not in entry]
        return self.organism_groups

    def resolve_organism_group(self, scientific_name: str) -> tuple[str, str]:
        scientific_name = scientific_name.strip()
        if not scientific_name or scientific_name.lower() == "unknown":
            return "", "No species assignment available for PD taxgroup selection."

        override = ORGANISM_GROUP_OVERRIDES.get(scientific_name.lower())
        if override:
            return override, ""

        groups = self.list_organism_groups()
        exact_name = scientific_name.replace(" ", "_")
        genus = scientific_name.split()[0]

        if exact_name in groups:
            return exact_name, ""

        exact_prefix_matches = [group for group in groups if group.startswith(f"{exact_name}_")]
        if len(exact_prefix_matches) == 1:
            return exact_prefix_matches[0], ""

        if genus in groups:
            return genus, ""

        genus_prefix_matches = [group for group in groups if group.startswith(f"{genus}_")]
        if len(genus_prefix_matches) == 1:
            return genus_prefix_matches[0], ""

        normalized_scientific = normalize_name(scientific_name)
        normalized_matches = [group for group in groups if normalize_name(group.replace("_", " ")) == normalized_scientific]
        if len(normalized_matches) == 1:
            return normalized_matches[0], ""

        return "", f"No PD organism group match found for '{scientific_name}'."

    def get_group_release(self, group: str) -> dict[str, str]:
        if group in self.group_releases:
            return self.group_releases[group]

        base_url = urljoin(FTP_RESULTS_URL, f"{group}/latest_snps/")
        metadata_entries = fetch_directory_entries(urljoin(base_url, "Metadata/"))
        cluster_entries = fetch_directory_entries(urljoin(base_url, "Clusters/"))
        exception_entries = fetch_directory_entries(urljoin(base_url, "Exceptions/"))

        metadata_name = next((entry for entry in metadata_entries if entry.endswith(".metadata.tsv")), "")
        all_isolates_name = next((entry for entry in cluster_entries if entry.endswith(".reference_target.all_isolates.tsv")), "")
        cluster_list_name = next((entry for entry in cluster_entries if entry.endswith(".reference_target.cluster_list.tsv")), "")
        snp_distances_name = next((entry for entry in cluster_entries if entry.endswith(".reference_target.SNP_distances.tsv")), "")
        exceptions_name = next((entry for entry in exception_entries if entry.endswith(".exceptions.tsv")), "")

        release = {
            "metadata_url": urljoin(base_url, f"Metadata/{metadata_name}") if metadata_name else "",
            "all_isolates_url": urljoin(base_url, f"Clusters/{all_isolates_name}") if all_isolates_name else "",
            "cluster_list_url": urljoin(base_url, f"Clusters/{cluster_list_name}") if cluster_list_name else "",
            "snp_distances_url": urljoin(base_url, f"Clusters/{snp_distances_name}") if snp_distances_name else "",
            "exceptions_url": urljoin(base_url, f"Exceptions/{exceptions_name}") if exceptions_name else "",
        }
        self.group_releases[group] = release
        return release


def select_query_rows(records: pd.DataFrame, release: dict[str, str]) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    sample_by_biosample = defaultdict(list)
    sample_by_run = defaultdict(list)
    sample_by_target = defaultdict(list)
    sample_by_asm = defaultdict(list)
    for row in records.itertuples(index=False):
        if row.biosample_acc:
            sample_by_biosample[row.biosample_acc].append(row.Sample)
        if row.srr_acc:
            sample_by_run[row.srr_acc].append(row.Sample)
        if row.pd_target_acc:
            sample_by_target[row.pd_target_acc].append(row.Sample)
        if row.pd_asm_acc:
            sample_by_asm[row.pd_asm_acc].append(row.Sample)

    query_rows: dict[str, dict[str, str]] = {}
    if release["metadata_url"]:
        for row in stream_tsv_rows(release["metadata_url"]):
            matches = []
            matches.extend(sample_by_biosample.get(row.get("biosample_acc", ""), []))
            matches.extend(sample_by_run.get(row.get("Run", ""), []))
            matches.extend(sample_by_target.get(row.get("target_acc", ""), []))
            matches.extend(sample_by_asm.get(row.get("asm_acc", ""), []))
            for sample in matches:
                query_rows.setdefault(sample, row)

    exception_rows: dict[str, dict[str, str]] = {}
    if release["exceptions_url"]:
        for row in stream_tsv_rows(release["exceptions_url"]):
            matches = []
            matches.extend(sample_by_biosample.get(row.get("biosample_acc", ""), []))
            for accession in re.split(r"[,\s;]+", row.get("run(s)", "") or row.get("run", "")):
                if accession:
                    matches.extend(sample_by_run.get(accession, []))
            for sample in matches:
                exception_rows.setdefault(sample, row)

    return query_rows, exception_rows


def load_cluster_assignments(targets: set[str], release: dict[str, str]) -> dict[str, dict[str, str]]:
    assignments = {}
    if not targets or not release["all_isolates_url"]:
        return assignments
    for row in stream_tsv_rows(release["all_isolates_url"]):
        target_acc = row.get("target_acc", "")
        if target_acc in targets:
            assignments[target_acc] = row
    return assignments


def load_cluster_members(clusters: set[str], release: dict[str, str]) -> dict[str, list[dict[str, str]]]:
    members = defaultdict(list)
    if not clusters or not release["cluster_list_url"]:
        return members
    for row in stream_tsv_rows(release["cluster_list_url"]):
        cluster = row.get("PDS_acc", "")
        if cluster in clusters:
            members[cluster].append(row)
    return members


def load_pairwise_distances(query_targets: set[str], clusters: set[str], release: dict[str, str]) -> tuple[dict[str, dict[str, Optional[int]]], str]:
    if not query_targets or not release["snp_distances_url"]:
        return {}, "No SNP distance table available for this taxgroup."

    size_bytes = get_remote_size_bytes(release["snp_distances_url"])
    if size_bytes is not None and size_bytes > FTP_SCANABLE_DISTANCE_BYTES:
        size_gb = round(size_bytes / (1024 ** 3), 1)
        return {}, f"SNP distance table is {size_gb} GB; using same-cluster comparators instead."

    distances: dict[str, dict[str, Optional[int]]] = defaultdict(dict)
    try:
        for row in stream_tsv_rows(release["snp_distances_url"]):
            cluster = row.get("PDS_acc", "")
            if clusters and cluster not in clusters:
                continue
            target_1 = row.get("target_acc_1", "")
            target_2 = row.get("target_acc_2", "")
            if target_1 in query_targets:
                other = target_2
                query_target = target_1
            elif target_2 in query_targets:
                other = target_1
                query_target = target_2
            else:
                continue

            distance = row.get("compatible_distance", "")
            try:
                numeric_distance = int(distance)
            except ValueError:
                numeric_distance = None

            previous = distances[query_target].get(other)
            if previous is None:
                distances[query_target][other] = numeric_distance
            elif numeric_distance is not None and (previous is None or numeric_distance < previous):
                distances[query_target][other] = numeric_distance
    except (HTTPError, URLError, TimeoutError) as exc:
        return {}, f"Unable to scan SNP distance table: {exc}"

    return distances, ""


def load_metadata_rows_by_target(targets: set[str], release: dict[str, str]) -> dict[str, dict[str, str]]:
    rows = {}
    if not targets or not release["metadata_url"]:
        return rows
    for row in stream_tsv_rows(release["metadata_url"]):
        target_acc = row.get("target_acc", "")
        if target_acc in targets:
            rows[target_acc] = row
    return rows


def note_from_exception(row: dict[str, str]) -> str:
    return row.get("exception", "") or row.get("consequence", "") or row.get("exception type", "")


def build_comparator_rows(
    sample: str,
    query_target: str,
    cluster: str,
    taxgroup: str,
    mode: str,
    selected_targets: list[str],
    distances: dict[str, Optional[int]],
    comparator_metadata: dict[str, dict[str, str]],
    cluster_members: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    for rank, target_acc in enumerate(selected_targets, start=1):
        metadata = comparator_metadata.get(target_acc, {})
        cluster_member = cluster_members.get(target_acc, {})
        distance = distances.get(target_acc)
        rows.append(
            {
                "Sample": sample,
                "pd_query_target_acc": query_target,
                "pd_snp_cluster": cluster,
                "pd_taxgroup": taxgroup,
                "comparator_rank": str(rank),
                "comparator_mode": mode,
                "comparator_distance": "" if distance is None else str(distance),
                "comparator_target_acc": target_acc,
                "comparator_biosample_acc": metadata.get("biosample_acc", "") or cluster_member.get("biosample_acc", ""),
                "comparator_asm_acc": metadata.get("asm_acc", "") or cluster_member.get("gencoll_acc", ""),
                "comparator_run": metadata.get("Run", ""),
                "comparator_sample_name": metadata.get("sample_name", ""),
                "comparator_scientific_name": metadata.get("scientific_name", ""),
                "comparator_source_type": metadata.get("source_type", ""),
                "comparator_host": metadata.get("host", ""),
                "comparator_isolation_source": metadata.get("isolation_source", ""),
            }
        )
    return rows


def run_live_lookup(records: pd.DataFrame, comparator_limit: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    client = PDLiveClient()
    comparator_rows: list[dict[str, str]] = []

    group_to_indices = defaultdict(list)
    for idx, row in records.iterrows():
        group, note = client.resolve_organism_group(row["scientific_name"])
        records.at[idx, "pd_taxgroup"] = group
        if group:
            group_to_indices[group].append(idx)
        elif row["biosample_acc"]:
            records.at[idx, "pd_lookup_status"] = "UNSUPPORTED_ORGANISM"
            records.at[idx, "pd_lookup_source"] = "pd_ftp_latest_snps"
            records.at[idx, "pd_lookup_note"] = note

    for group, indices in group_to_indices.items():
        release = client.get_group_release(group)
        group_records = records.loc[indices, ["Sample", "srr_acc", "biosample_acc", "pd_target_acc", "pd_asm_acc"]]
        query_rows, exception_rows = select_query_rows(group_records, release)

        query_targets = {
            row.get("target_acc", "") for row in query_rows.values() if row.get("target_acc", "")
        } | {
            records.at[idx, "pd_target_acc"] for idx in indices if records.at[idx, "pd_target_acc"]
        }
        assignments = load_cluster_assignments(query_targets, release)
        clusters = {row.get("PDS_acc", "") for row in assignments.values() if row.get("PDS_acc", "")}
        cluster_members = load_cluster_members(clusters, release)
        pairwise_distances, pairwise_note = load_pairwise_distances(set(assignments.keys()), clusters, release)

        comparator_target_lookup: dict[str, list[str]] = {}
        comparator_mode_lookup: dict[str, str] = {}
        comparator_distance_lookup: dict[str, dict[str, Optional[int]]] = {}
        needed_comparator_targets: set[str] = set()

        for idx in indices:
            sample = records.at[idx, "Sample"]
            query_row = query_rows.get(sample)
            exception_row = exception_rows.get(sample)

            if query_row:
                records.at[idx, "pd_target_acc"] = query_row.get("target_acc", "") or records.at[idx, "pd_target_acc"]
                records.at[idx, "pd_asm_acc"] = query_row.get("asm_acc", "") or records.at[idx, "pd_asm_acc"]
                records.at[idx, "biosample_acc"] = query_row.get("biosample_acc", "") or records.at[idx, "biosample_acc"]
                records.at[idx, "pd_scientific_name"] = query_row.get("scientific_name", "")
                records.at[idx, "pd_source_type"] = query_row.get("source_type", "")
                records.at[idx, "pd_host"] = query_row.get("host", "")
                records.at[idx, "pd_isolation_source"] = query_row.get("isolation_source", "")
                records.at[idx, "pd_lookup_source"] = "pd_ftp_latest_snps"

            target_acc = records.at[idx, "pd_target_acc"]
            assignment = assignments.get(target_acc, {})
            cluster = assignment.get("PDS_acc", "")
            records.at[idx, "pd_snp_cluster"] = cluster

            if query_row and cluster:
                distance_hits = pairwise_distances.get(target_acc, {})
                if distance_hits:
                    selected_targets = [
                        other_target
                        for other_target, _distance in sorted(
                            distance_hits.items(),
                            key=lambda item: (item[1] is None, item[1] if item[1] is not None else 10**9, item[0]),
                        )
                        if other_target != target_acc
                    ][:comparator_limit]
                    comparator_mode = "pairwise_distance"
                    comparator_note = pairwise_note or "Comparators ranked by PD compatible SNP distance."
                else:
                    selected_targets = [
                        member.get("target_acc", "")
                        for member in cluster_members.get(cluster, [])
                        if member.get("target_acc", "") and member.get("target_acc", "") != target_acc
                    ][:comparator_limit]
                    comparator_mode = "same_cluster"
                    comparator_note = pairwise_note or "Comparators selected from the same current PD SNP cluster."

                comparator_target_lookup[sample] = selected_targets
                comparator_mode_lookup[sample] = comparator_mode
                comparator_distance_lookup[sample] = distance_hits
                needed_comparator_targets.update(selected_targets)

                records.at[idx, "pd_comparator_mode"] = comparator_mode
                records.at[idx, "pd_comparator_count"] = str(len(selected_targets))
                records.at[idx, "pd_lookup_status"] = "FOUND"
                records.at[idx, "pd_lookup_note"] = comparator_note
            elif query_row and not cluster:
                records.at[idx, "pd_lookup_status"] = "FOUND_NO_CLUSTER"
                records.at[idx, "pd_lookup_note"] = "PD metadata record found, but no current SNP-cluster assignment was published."
            elif exception_row:
                records.at[idx, "pd_lookup_status"] = "QC_EXCEPTION"
                records.at[idx, "pd_lookup_source"] = "pd_ftp_exceptions"
                records.at[idx, "pd_lookup_note"] = note_from_exception(exception_row)
            elif records.at[idx, "biosample_acc"]:
                records.at[idx, "pd_lookup_status"] = "NOT_FOUND"
                records.at[idx, "pd_lookup_source"] = "pd_ftp_latest_snps"
                records.at[idx, "pd_lookup_note"] = "BioSample was not present in the current public PD FTP release."

        comparator_metadata = load_metadata_rows_by_target(needed_comparator_targets, release)
        member_lookup = {
            cluster: {member.get("target_acc", ""): member for member in members if member.get("target_acc", "")}
            for cluster, members in cluster_members.items()
        }

        for idx in indices:
            sample = records.at[idx, "Sample"]
            selected_targets = comparator_target_lookup.get(sample, [])
            if not selected_targets:
                continue
            target_acc = records.at[idx, "pd_target_acc"]
            cluster = records.at[idx, "pd_snp_cluster"]
            comparator_rows.extend(
                build_comparator_rows(
                    sample=sample,
                    query_target=target_acc,
                    cluster=cluster,
                    taxgroup=group,
                    mode=comparator_mode_lookup.get(sample, "same_cluster"),
                    selected_targets=selected_targets,
                    distances=comparator_distance_lookup.get(sample, {}),
                    comparator_metadata=comparator_metadata,
                    cluster_members=member_lookup.get(cluster, {}),
                )
            )

    return records[WORKING_COLUMNS], pd.DataFrame(comparator_rows, columns=COMPARATOR_COLUMNS).fillna("")


def run_table_lookup(
    records: pd.DataFrame,
    pd_isolates_path: Optional[Path],
    pd_exceptions_path: Optional[Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pd_isolates, isolates_error = load_pd_table(pd_isolates_path, PD_ISOLATE_ALIASES)
    pd_exceptions, exceptions_error = load_pd_table(pd_exceptions_path, PD_EXCEPTION_ALIASES)

    if pd_isolates is None:
        records["pd_lookup_status"] = "CONFIG_ERROR"
        records["pd_lookup_source"] = "pd_isolates_tsv"
        records["pd_lookup_note"] = isolates_error or "pd_isolates_tsv is required when pd_backend=table"
        return records[WORKING_COLUMNS], pd.DataFrame(columns=COMPARATOR_COLUMNS)

    found = records.merge(pd_isolates, on="biosample_acc", how="left", suffixes=("", "_pd")).fillna("")
    for column in [
        "pd_target_acc",
        "pd_asm_acc",
        "pd_snp_cluster",
        "pd_taxgroup",
        "pd_scientific_name",
        "pd_source_type",
        "pd_host",
        "pd_isolation_source",
    ]:
        pd_column = f"{column}_pd"
        if pd_column in found.columns:
            found[column] = found[pd_column].where(found[pd_column] != "", found[column])
    records = found[WORKING_COLUMNS].copy()

    if pd_exceptions is not None:
        exception_notes = pd_exceptions.rename(columns={"pd_lookup_note": "pd_exception_note"})
        records = records.merge(exception_notes[["biosample_acc", "srr_acc", "pd_exception_note"]], on=["biosample_acc", "srr_acc"], how="left").fillna("")
    else:
        records["pd_exception_note"] = ""

    if exceptions_error:
        records["pd_lookup_note"] = records["pd_lookup_note"].mask(records["pd_lookup_note"] == "", exceptions_error)

    def classify(row: pd.Series) -> tuple[str, str, str]:
        if row["pd_snp_cluster"]:
            return "FOUND", "pd_isolates_tsv", row["pd_lookup_note"]
        if row["pd_exception_note"]:
            return "QC_EXCEPTION", "pd_exceptions_tsv", row["pd_exception_note"]
        if not row["biosample_acc"]:
            if row["srr_acc"]:
                return "LOOKUP_ERROR", row["pd_lookup_source"] or "sra_run_new", row["pd_lookup_note"]
            return "NO_ACCESSION", "", "Provide BioSample/SRR in pd_sample_metadata_tsv or use SRR-like sample names."
        return "NOT_FOUND", "pd_isolates_tsv", row["pd_lookup_note"] or "BioSample not found in provided PD table."

    statuses = records.apply(classify, axis=1, result_type="expand")
    statuses.columns = ["pd_lookup_status", "pd_lookup_source", "pd_lookup_note"]
    records[["pd_lookup_status", "pd_lookup_source", "pd_lookup_note"]] = statuses
    return records[WORKING_COLUMNS], pd.DataFrame(columns=COMPARATOR_COLUMNS)


samples = [str(sample) for sample in snakemake.params["samples"]]
enabled = bool(snakemake.params["enabled"])
backend = str(snakemake.params.get("backend", "ftp")).strip().lower()
comparator_limit = int(snakemake.params.get("comparator_limit", 10))
sample_metadata_path = as_optional_path(snakemake.params["sample_metadata_tsv"])
pd_isolates_path = as_optional_path(snakemake.params["pd_isolates_tsv"])
pd_exceptions_path = as_optional_path(snakemake.params["pd_exceptions_tsv"])
mash_path = Path(str(snakemake.input["mash"]))
metadata_output_path = Path(str(snakemake.output["metadata"]))
comparator_output_path = Path(str(snakemake.output["comparators"]))
metadata_output_path.parent.mkdir(parents=True, exist_ok=True)

records = build_sample_records(samples, sample_metadata_path, mash_path)
for column in WORKING_COLUMNS:
    if column not in records.columns:
        records[column] = ""
records = records[WORKING_COLUMNS].copy()

comparator_records = pd.DataFrame(columns=COMPARATOR_COLUMNS)

if not enabled:
    records["pd_lookup_status"] = "DISABLED"
    records["pd_lookup_note"] = "Set pd_lookup=true to enable optional Pathogen Detection enrichment."
else:
    sra_cache: dict[str, dict[str, str]] = {}
    for idx, row in records.iterrows():
        if row["biosample_acc"] and row["scientific_name"]:
            continue
        if not row["srr_acc"]:
            continue
        sra_metadata = sra_cache.get(row["srr_acc"])
        if sra_metadata is None:
            sra_metadata = fetch_sra_metadata(row["srr_acc"])
            sra_cache[row["srr_acc"]] = sra_metadata
        if not row["biosample_acc"] and sra_metadata["biosample_acc"]:
            records.at[idx, "biosample_acc"] = sra_metadata["biosample_acc"]
        if not row.get("pd_lookup_note") and sra_metadata["note"]:
            records.at[idx, "pd_lookup_note"] = sra_metadata["note"]
        if not records.at[idx, "scientific_name"]:
            records.at[idx, "scientific_name"] = sra_metadata["scientific_name"]
        if sra_metadata["biosample_acc"]:
            records.at[idx, "pd_lookup_source"] = "sra_run_new"

    if backend in {"ftp", "live", "online"}:
        records, comparator_records = run_live_lookup(records, comparator_limit)
    elif backend in {"table", "local", "export"}:
        records, comparator_records = run_table_lookup(records, pd_isolates_path, pd_exceptions_path)
    else:
        records["pd_lookup_status"] = "CONFIG_ERROR"
        records["pd_lookup_source"] = "pd_backend"
        records["pd_lookup_note"] = f"Unsupported pd_backend '{backend}'. Use 'ftp' or 'table'."

    for idx, row in records.iterrows():
        if row["pd_lookup_status"]:
            continue
        if not row["biosample_acc"]:
            if row["srr_acc"]:
                records.at[idx, "pd_lookup_status"] = "LOOKUP_ERROR"
                records.at[idx, "pd_lookup_source"] = row["pd_lookup_source"] or "sra_run_new"
                if not row["pd_lookup_note"]:
                    records.at[idx, "pd_lookup_note"] = "Unable to resolve this run to a BioSample accession."
            else:
                records.at[idx, "pd_lookup_status"] = "NO_ACCESSION"
                records.at[idx, "pd_lookup_note"] = "Provide BioSample/SRR in pd_sample_metadata_tsv or use SRR-like sample names."
        else:
            records.at[idx, "pd_lookup_status"] = "NOT_FOUND"
            records.at[idx, "pd_lookup_source"] = row["pd_lookup_source"] or f"pd_{backend}"
            if not row["pd_lookup_note"]:
                records.at[idx, "pd_lookup_note"] = "No current PD record was found for this accession."

records = records[OUTPUT_COLUMNS].fillna("")
for column in COMPARATOR_COLUMNS:
    if column not in comparator_records.columns:
        comparator_records[column] = ""
comparator_records = comparator_records[COMPARATOR_COLUMNS].fillna("")

records.to_csv(metadata_output_path, sep="\t", index=False)
comparator_records.to_csv(comparator_output_path, sep="\t", index=False)
