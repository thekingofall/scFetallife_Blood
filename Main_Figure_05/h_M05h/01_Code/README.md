# Figure 5h: DNTT expression in liver pro-B cells

## Scripts

### `01_plot_dntt_expression_by_age.py`

Selects thymic DN(Q) T cells and liver pro-B cells from the atlas, groups them by developmental age, and plots `DNTT` expression. Dot color represents mean expression and dot size represents the fraction of expressing cells. It saves the Figure 4h/5h dot plots and selected-cell metadata.

## Data

**`M05h_DNTT_dotplot_data.tsv`**

Cell counts, mean DNTT expression, and the fraction of expressing cells for each age/population group.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
