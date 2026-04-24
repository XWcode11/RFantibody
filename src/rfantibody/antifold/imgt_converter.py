"""
HLT to AntiFold input converter.

Converts RFantibody HLT-format PDB files to AntiFold-compatible input:
- Removes target (T) chain
- Keeps only H and L chains
- Renumbers residues starting from 1 per chain (basic IMGT-like numbering)
- Maps CDR regions from HLT REMARKs to AntiFold IMGT region names
"""


class HLTtoAntiFoldConverter:
    """Convert HLT-format PDB to AntiFold input format."""

    def convert(self, pdb_path: str, heavy_chain: str = "H", light_chain: str = "L") -> list[str]:
        """Convert HLT PDB to AntiFold-compatible PDB lines.

        Args:
            pdb_path: Path to HLT-format PDB file
            heavy_chain: Chain ID for heavy chain (default "H")
            light_chain: Chain ID for light chain (default "L")

        Returns:
            List of PDB lines (H and L chains only, renumbered)
        """
        allowed_chains = {heavy_chain, light_chain}
        output_lines = []

        with open(pdb_path, "r") as f:
            lines = f.readlines()

        # First pass: collect residue numbers per chain to build renumbering maps
        chain_residues = {heavy_chain: [], light_chain: []}
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21:22].strip()
                if chain in allowed_chains:
                    resno = int(line[22:26].strip())
                    if resno not in chain_residues[chain]:
                        chain_residues[chain].append(resno)

        # Build renumbering maps: old_resno -> new_resno (1-indexed)
        renumber_map = {}
        for chain in allowed_chains:
            sorted_res = sorted(chain_residues[chain])
            renumber_map[chain] = {old: i + 1 for i, old in enumerate(sorted_res)}

        # Second pass: write filtered and renumbered lines
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21:22].strip()
                if chain not in allowed_chains:
                    continue
                # Renumber residue
                old_resno = int(line[22:26].strip())
                new_resno = renumber_map[chain].get(old_resno, old_resno)
                new_line = line[:22] + f"{new_resno:4d}" + line[26:]
                output_lines.append(new_line)
            elif line.startswith("TER"):
                chain = line[21:22].strip() if len(line) > 21 else ""
                if chain in allowed_chains or chain == "":
                    output_lines.append(line)
            elif line.startswith("REMARK PDBinfo-LABEL"):
                # Only keep CDR remarks for H or L chains
                # REMARK format: REMARK PDBinfo-LABEL:  <resi> <loop>
                parts = line.strip().split()
                if len(parts) >= 4:
                    loop_name = parts[-1]
                    if loop_name.startswith("H") or loop_name.startswith("L"):
                        output_lines.append(line)
            else:
                # Keep other non-ATOM lines (HEADER, TITLE, etc.)
                output_lines.append(line)

        return output_lines

    def infer_antifold_regions(self, pdb_path: str) -> list[str]:
        """Infer AntiFold design regions from HLT REMARK PDBinfo-LABEL lines.

        Args:
            pdb_path: Path to HLT-format PDB file

        Returns:
            List of AntiFold region names (e.g., ["CDRH1", "CDRH2", "CDRH3"])
        """
        hlt_loops = set()
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("REMARK PDBinfo-LABEL"):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        loop_name = parts[-1].upper()
                        if loop_name in ("H1", "H2", "H3", "L1", "L2", "L3"):
                            hlt_loops.add(loop_name)

        # Map HLT loop names to AntiFold IMGT region names
        mapping = {
            "H1": "CDRH1",
            "H2": "CDRH2",
            "H3": "CDRH3",
            "L1": "CDRL1",
            "L2": "CDRL2",
            "L3": "CDRL3",
        }

        regions = [mapping[loop] for loop in sorted(hlt_loops) if loop in mapping]
        return regions
