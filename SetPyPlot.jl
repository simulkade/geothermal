module SetPyPlot

# setting pyplot parameters to default values
using PyPlot

export setrcparam, plotyy

"""
setrcparam(;tick_font_size = 8,
legend_font_size = 8,
axis_label_font_size = 9,
fig_width = 9,
fig_height = 7)

sets the font size parameters to my favorite default values
"""
function setrcparam(;tick_font_size = 8,
    legend_font_size = 8,
    axis_label_font_size = 9,
    fig_width = 9,
    fig_height = 7)

    rc("xtick", labelsize=tick_font_size)
    rc("ytick", labelsize=tick_font_size)
    rc("legend", fontsize=legend_font_size)
    rc("axes", labelsize=axis_label_font_size)
    rc("figure", dpi = 300)
    rc("figure", figsize = (fig_width, fig_height)./2.54)
    rc("savefig", bbox = "tight")
end

"""
function plotyy(x1, y1, x2, y2; x_label = "", y1_label = "", y2_label = "", y1_style = "-bo", y2_style = "--vr")
"""
function plotyy(x1, y1, x2, y2; fig_size = (3.5, 2.5), x_label = "", y1_label = "", y2_label = "", y1_style = "-bo", y2_style = "--vr")
    fig, ax1 = subplots(figsize = fig_size)
    ax1.plot(x1, y1, y1_style)
    ax1.set_ylabel(y1_label, color="b")
    ax1.set_xlabel(x_label, color="b")
    ax1.tick_params("y", colors = "b")

    ax2 = ax1.twinx()
    ax2.plot(x2, y2, y2_style)
    ax2.set_ylabel(y2_label, color="r")
    ax2.tick_params("y", colors = "r")

    tight_layout()

    return fig, ax1, ax2
end
# more to come later...

end # module
