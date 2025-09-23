# Functions dispaying metadata in plots
import matplotlib.pyplot as plt

from src.color_maps import all_color_maps

plt.rcParams.update({"font.family": "DejaVu Sans"})
PLT_STYLE = "tableau-colorblind10"
plt.style.use(PLT_STYLE)


def display_diet_information(df, diet_var, age_var, row_label, title="", fonts=14):
    if diet_var == "diet_weaning":
        group_order = ["no", "yes", "unknown"]
    elif diet_var == "diet_milk":
        group_order = ["bd", "mixed", "fd", "unknown"]

    # Copy and enforce "unknown" color
    color_dict = {**all_color_maps[diet_var]}
    color_dict["unknown"] = "#d9d9d8"

    test_df = df[[age_var, diet_var]].copy()
    test_df[diet_var] = test_df[diet_var].astype(str)
    test_df[diet_var] = test_df[diet_var].replace("nan", "unknown")

    df_grouped = test_df.groupby([age_var, diet_var]).size().unstack(fill_value=0)
    df_grouped = df_grouped.sort_index()
    df_grouped = df_grouped.reindex(columns=group_order)

    ax = df_grouped.plot.bar(
        stacked=True,
        figsize=(9, 2.5),
        color=[color_dict.get(col, "#d9d9d8") for col in df_grouped.columns],
    )

    # Only show x-axis labels for specific ages, horizontally
    desired_ticks = [5, 10, 15, 20, 25, 30, 35]
    index_vals = list(df_grouped.index)
    pos_map = {val: i for i, val in enumerate(index_vals)}
    tick_positions = [pos_map[v] for v in desired_ticks if v in pos_map]
    tick_labels = [str(v) for v in desired_ticks if v in pos_map]
    if tick_positions:
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=0, ha="center")
    else:
        # Ensure horizontal labels even if nothing matched
        for label in ax.get_xticklabels():
            label.set_rotation(0)

    ax.set_xlabel("Age [months]", fontsize=fonts)
    ax.set_ylabel(f"# {row_label}", fontsize=fonts)
    ax.set_title(title, fontsize=fonts)
    ax.tick_params(axis="both", labelsize=fonts)
    legend = ax.get_legend()
    if legend is not None:
        legend.set_title(legend.get_title().get_text(), prop={"size": fonts})
        for txt in legend.get_texts():
            txt.set_fontsize(fonts)

    return ax
