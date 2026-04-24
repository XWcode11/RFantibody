"""
HLT to AntiFold input converter.

Converts RFantibody HLT-format PDB files to AntiFold-compatible input:
- Removes target (T) chain, keeps only H and L chains
- Renumbers residues so CDR positions align to AntiFold's IMGT anchors
  (CDR1 at 27, CDR2 at 56, CDR3 at 105); framework residues fill in from
  CDR anchors outward, matching the IMGT_dict ranges used by AntiFold
- Maps CDR regions from HLT REMARKs to AntiFold IMGT region names
"""

_CDR_ANCHORS = {"1": 27, "2": 56, "3": 105}
_FW4_ANCHOR = 118
_FW1_END = 26


class HLTtoAntiFoldConverter:
    """Convert HLT-format PDB to AntiFold input format."""

    def convert(
        self,
        pdb_path: str,
        heavy_chain: str = "H",
        light_chain: str = "L",
    ) -> list[str]:
        """Convert HLT PDB to AntiFold-compatible PDB lines."""
        allowed_chains = {heavy_chain, light_chain}

        with open(pdb_path, "r") as f:
            lines = f.readlines()

        chain_resnos: dict[str, list[int]] = {heavy_chain: [], light_chain: []}
        cdr_labels: dict[str, dict[int, str]] = {heavy_chain: {}, light_chain: {}}

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21:22].strip()
                if chain in allowed_chains:
                    resno = int(line[22:26].strip())
                    if resno not in chain_resnos[chain]:
                        chain_resnos[chain].append(resno)
            elif line.startswith("REMARK PDBinfo-LABEL"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        resno = int(parts[2])
                    except ValueError:
                        continue
                    label = parts[3]
                    if label[:1] in ("H", "L") and label[1:] in ("1", "2", "3"):
                        prefix = label[0]
                        if prefix in cdr_labels:
                            cdr_labels[prefix][resno] = label

        renumber_map: dict[str, dict[int, int]] = {}
        for chain in allowed_chains:
            renumber_map[chain] = self._build_imgt_numbering(
                sorted(chain_resnos[chain]),
                cdr_labels[chain],
            )

        output_lines: list[str] = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21:22].strip()
                if chain not in allowed_chains:
                    continue
                old_resno = int(line[22:26].strip())
                new_resno = renumber_map[chain].get(old_resno, old_resno)
                output_lines.append(line[:22] + f"{new_resno:4d}" + line[26:])
            elif line.startswith("TER"):
                chain = line[21:22].strip() if len(line) > 21 else ""
                if chain in allowed_chains or chain == "":
                    output_lines.append(line)
            elif line.startswith("REMARK PDBinfo-LABEL"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    loop_name = parts[-1]
                    if loop_name.startswith("H") or loop_name.startswith("L"):
                        output_lines.append(line)
            else:
                output_lines.append(line)

        return output_lines

    @staticmethod
    def _build_imgt_numbering(
        resnos: list[int],
        cdr_labels: dict[int, str],
    ) -> dict[int, int]:
        """Build old->new residue map aligning CDRs to IMGT anchors 27/56/105."""
        if not resnos:
            return {}

        cdrs: dict[str, list[int]] = {"1": [], "2": [], "3": []}
        for r in sorted(cdr_labels.keys()):
            key = cdr_labels[r][1]
            if key in cdrs:
                cdrs[key].append(r)

        mapping: dict[int, int] = {}
        cdr_new_end: dict[str, int | None] = {"1": None, "2": None, "3": None}
        for k in ("1", "2", "3"):
            for i, r in enumerate(cdrs[k]):
                mapping[r] = _CDR_ANCHORS[k] + i
            if cdrs[k]:
                cdr_new_end[k] = _CDR_ANCHORS[k] + len(cdrs[k]) - 1

        first_cdr_old = None
        first_cdr_key = None
        for k in ("1", "2", "3"):
            if cdrs[k]:
                first_cdr_old = cdrs[k][0]
                first_cdr_key = k
                break

        if first_cdr_old is None:
            for i, r in enumerate(resnos):
                mapping[r] = i + 1
            return mapping

        fw1 = [r for r in resnos if r < first_cdr_old]
        fw1_target_end = _FW1_END if cdrs["1"] else _CDR_ANCHORS[first_cdr_key] - 1
        start = max(1, fw1_target_end - len(fw1) + 1)
        for i, r in enumerate(fw1):
            mapping[r] = start + i

        def _fill_middle(prev_key: str, next_key: str) -> None:
            if not cdrs[prev_key] or not cdrs[next_key]:
                return
            lo = cdrs[prev_key][-1]
            hi = cdrs[next_key][0]
            block = [r for r in resnos if lo < r < hi]
            prev_end = cdr_new_end[prev_key]
            next_start = _CDR_ANCHORS[next_key]
            packed_start = next_start - len(block)
            if packed_start <= prev_end:
                start_i = prev_end + 1
                for i, r in enumerate(block):
                    mapping[r] = start_i + i
            else:
                for i, r in enumerate(block):
                    mapping[r] = packed_start + i

        _fill_middle("1", "2")
        _fill_middle("2", "3")

        last_cdr_key = next(k for k in ("3", "2", "1") if cdrs[k])
        last_cdr_old = cdrs[last_cdr_key][-1]
        fw4 = [r for r in resnos if r > last_cdr_old]
        prev_end = cdr_new_end[last_cdr_key]
        start_fw4 = max(_FW4_ANCHOR, prev_end + 1) if last_cdr_key == "3" else prev_end + 1
        for i, r in enumerate(fw4):
            mapping[r] = start_fw4 + i

        return mapping

    def infer_antifold_regions(self, pdb_path: str) -> list[str]:
        """Infer AntiFold design regions from HLT REMARK PDBinfo-LABEL lines."""
        hlt_loops: set[str] = set()
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("REMARK PDBinfo-LABEL"):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        loop_name = parts[-1].upper()
                        if loop_name in ("H1", "H2", "H3", "L1", "L2", "L3"):
                            hlt_loops.add(loop_name)

        mapping = {
            "H1": "CDRH1", "H2": "CDRH2", "H3": "CDRH3",
            "L1": "CDRL1", "L2": "CDRL2", "L3": "CDRL3",
        }
        return [mapping[loop] for loop in sorted(hlt_loops) if loop in mapping]
