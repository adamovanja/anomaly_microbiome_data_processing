# Functions dispaying metadata in plots
import matplotlib.pyplot as plt
import seaborn as sns

from src.color_maps import all_color_maps

plt.rcParams.update({"font.family": "DejaVu Sans"})
plt.rcParams.update({"font.size": 10})
PLT_STYLE = "tableau-colorblind10"
plt.style.use(PLT_STYLE)


def display_diet_information(df, diet_var, age_var, row_label, title=""):
    if diet_var == "diet_weaning":
        group_order = ["no", "yes", "finished", "unknown"]
    elif diet_var == "diet_milk":
        group_order = ["bd", "mixed", "fd", "no milk", "unknown"]
    color_dict = all_color_maps[diet_var]
    test_df = df[[age_var, diet_var]].copy()
    test_df[diet_var] = test_df[diet_var].astype(str).replace("nan", "unknown")

    df_grouped = test_df.groupby([age_var, diet_var]).size().unstack(fill_value=0)
    df_grouped = df_grouped.sort_index().reindex(columns=group_order)
    ax = df_grouped.plot.bar(
        stacked=True,
        figsize=(20, 5),
        color=[color_dict[col] for col in df_grouped.columns],
    )
    ax.set_xlabel("Age [months]" if age_var == "age_months_rounded1" else age_var)
    ax.set_xticklabels([str(int(x)) for x in df_grouped.index], rotation=0, ha="center")
    ax.set_ylabel(f"# {row_label}")
    ax.set_title(title, fontsize=10, y=1.0)
    return ax


def display_diet_information_in_one(df, diet_vars, age_var, row_label, title=""):
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 3), dpi=400)

    for i, dv in enumerate(diet_vars):
        if dv == "diet_weaning":
            group_order = ["no", "yes", "unknown"]
        else:
            group_order = ["bd", "mixed", "fd", "unknown"]
        color_dict = all_color_maps[dv]
        test_df = df[[age_var, dv]].copy()
        test_df[dv] = test_df[dv].astype(str).replace("nan", "unknown")

        df_grouped = test_df.groupby([age_var, dv]).size().unstack(fill_value=0)
        df_grouped = df_grouped.sort_index().reindex(columns=group_order)

        ax = df_grouped.plot.bar(
            ax=axes[i],
            stacked=True,
            color=[color_dict[col] for col in df_grouped.columns],
        )
        ax.set_xlabel("Age [months]" if age_var == "age_months_rounded1" else age_var)
        ax.set_xticklabels(
            [str(int(x)) for x in df_grouped.index], rotation=0, ha="center", fontsize=8
        )
        ax.set_ylabel(f"# {row_label}")
        handles, labels = ax.get_legend_handles_labels()
        legend_title = ax.get_legend().get_title().get_text() if ax.get_legend() else ""

        ax.legend(
            handles,
            labels,
            title=legend_title.replace("_", " "),
            loc="upper right",
            bbox_to_anchor=(1.18, 1.12),
        )
    fig.suptitle(title, fontsize=12, y=0.9)
    fig.tight_layout()
    return fig, axes


def plot_box_violin(y, color, ax, horizontal=False):
    sns.violinplot(
        x=y if horizontal else None,
        y=y if not horizontal else None,
        inner=None,
        ax=ax,
        color=color,
        orient="h" if horizontal else "v",
    )
    sns.boxplot(
        x=y if horizontal else None,
        y=y if not horizontal else None,
        width=0.1,
        boxprops={"facecolor": "white", "edgecolor": "black", "zorder": 2},
        flierprops={
            "marker": "o",
            "markerfacecolor": "none",
            "markeredgecolor": "black",
        },
        ax=ax,
        orient="h" if horizontal else "v",
    )
