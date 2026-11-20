import gzip

from entropy_age.experiments.bj_telomerase import HTERT, MORTAL_EARLIER, MORTAL_SAME, groups, read_cultures


def _matrix(tmp_path, rows):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("!Sample_title\t" + "\t".join(f'"{r[0]}"' for r in rows) + "\n")
        handle.write("!Sample_geo_accession\t" + "\t".join(f'"GSM{i}"' for i in range(len(rows))) + "\n")
        handle.write(
            "!Sample_supplementary_file\t"
            + "\t".join(f'"ftp://x/GSM{i}_{r[1]}_R01C01_Grn.idat.gz"' for i, r in enumerate(rows))
            + "\n"
        )
        handle.write(
            "!Sample_supplementary_file\t"
            + "\t".join(f'"ftp://x/GSM{i}_{r[1]}_R01C01_Red.idat.gz"' for i, r in enumerate(rows))
            + "\n"
        )
    return path


def test_read_cultures_and_groups(tmp_path):
    rows = [
        ("Cell_BJF_hTERT_Trajectory_B_PDL_38", "207222490001"),
        ("Cell_BJF_Trajectory_A_PDL_22", "207222490001"),
        ("Cell_BJF_Trajectory_A_PDL_2.3_Replicate_1", "206439850165"),
        ("Cell_BJF_Trajectory_DeepSen_PDL_104_Replicate_1", "206439850165"),
    ]
    cultures = read_cultures(_matrix(tmp_path, rows))
    assert cultures["telomerase"].to_list() == [True, False, False, False]
    assert cultures["doublings"].to_list() == [38.0, 22.0, 2.3, 104.0]
    assert cultures["chip"].to_list() == ["207222490001", "207222490001", "206439850165", "206439850165"]
    g = groups(cultures)
    assert g[HTERT]["sample"].to_list() == ["GSM0"]
    assert g[MORTAL_SAME]["sample"].to_list() == ["GSM1"]
    assert g[MORTAL_EARLIER]["sample"].to_list() == ["GSM2"]
