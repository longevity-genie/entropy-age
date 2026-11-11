"""Reproductions of the paper's figures and the extension experiments."""

from collections.abc import Callable
from pathlib import Path

from entropy_age.experiments import (
    atlas_real_data,
    bxd_liver,
    chromatin_clocks,
    chromatin_match,
    clock_cpgs,
    consortium_interventions,
    extension_selection,
    fig1_bulk,
    fig3_single_cell,
    fig4_published_clocks,
    fig5_real_data,
    gene_sets,
    human_tissues,
    intergenic,
    interventions,
    mouse_rrbs_tissues,
    mouse_tissues,
    passaged_cells,
    per_person,
    rna_gtex,
    rna_rodents,
    species_shield,
    tissue_contrast,
)
from entropy_age.experiments.common import Context, Scale, ScaleName

FIGURES: dict[str, Callable[[Context], list[Path]]] = {
    "fig1": fig1_bulk.run,
    "fig3": fig3_single_cell.run,
    "fig4": fig4_published_clocks.run,
    "fig5": fig5_real_data.run,
    "extension": extension_selection.run,
    "atlas": atlas_real_data.run,
    "genesets": gene_sets.run,
    "tissues": tissue_contrast.run,
    "tissues_cells": tissue_contrast.run_cells,
    "mouse": mouse_tissues.run,
    "gtex": human_tissues.run,
    "fibres": human_tissues.run_fibres,
    "rna_gtex": rna_gtex.run,
    "rna_rodents": rna_rodents.run,
    "passage": passaged_cells.run,
    "passage_sturm": passaged_cells.run_sturm,
    "chromatin": chromatin_match.run,
    "chromatin_clocks": chromatin_clocks.run,
    "interventions": interventions.run,
    "bxd": bxd_liver.run,
    "rrbs_tissues": mouse_rrbs_tissues.run,
    "species": species_shield.run,
    "mmc_interventions": consortium_interventions.run,
    "intergenic": intergenic.run,
    "clocks": clock_cpgs.run,
    "per_person": per_person.run,
}

__all__ = ["FIGURES", "Context", "Scale", "ScaleName"]
