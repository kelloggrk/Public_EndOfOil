# ==============================================================
# This file holds functions for plotting results from the model
# ==============================================================



# ==============================================================
# Function for stacked drilling, production, and price plots
# ==============================================================
function comboplot(Oi::NamedTuple, Od::NamedTuple, Om::NamedTuple, 
    T::Int64, serif::Bool = false, L::Int64 = 3)
    # INPUTS
    # Oi: tuple of baseline outcomes
    # Od: tuple of anticipated demand decline outcomes
    # Om: tuple of myopic demand decline outcomes
    # T: number of periods per year
    # serif: true/false flag for whether to use a serif font
    # L: number of lines to plot (baseline, anticipated decline, unanticipated decline)

    # OUTPUTS
    # combined_plot: plot of drilling, production, and price

    # First set up time axis and region labels
    FirstY = 2023; LastY = 2119;        # first and last year to plot
    NTplot = (LastY-FirstY) * T + 1;  # number of times to plot
    TimeY = [FirstY:1/T:LastY];       # vector of time in calendar years
 
    # Extract aggregate drilling, production, and price data
    # Set drilling to zero for models with no investment
    if haskey(Oi, :dvec) == false
        Y11 = zeros(NTplot)
        Y12 = zeros(NTplot)
        Y13 = zeros(NTplot)
    else
        Y11 = vec(sum(Oi.dvec[1:NTplot,:]*T, dims=2));    # total drilling in all regions
        Y12 = vec(sum(Od.dvec[1:NTplot,:]*T, dims=2));
        Y13 = vec(sum(Om.dvec[1:NTplot,:]*T, dims=2));  
    end
    Y21 = vec(sum(Oi.qvec[1:NTplot,:], dims=2));        # total prod in all regions
    Y22 = vec(sum(Od.qvec[1:NTplot,:], dims=2)); 
    Y23 = vec(sum(Om.qvec[1:NTplot,:], dims=2));  
    Y31 = vec(sum(Oi.pvec[1:NTplot,:], dims=2));        # oil price
    Y32 = vec(sum(Od.pvec[1:NTplot,:], dims=2));  
    Y33 = vec(sum(Om.pvec[1:NTplot,:], dims=2)); 
    
    # More plot variables
    Tgap = 20;                          # gap between x-axis labels
    maxD = max(10,round(maximum(Y11)/10)*10); maxQ = 100; maxP = 251   # max Y axis values
    fsize = 20;                         # font size 
    fsizess = 23;                       # font size for slides

    
    # Legend labels
    strb = "baseline demand"
    stra = "anticipated demand decline"
    stru = "unanticipated demand decline"

    # Plot drilling
    dplot = plot(TimeY, Y11, label=strb, linewidth=1.5, linecolor=:red, linestyle=:solid)
    if L>=2
        plot!(TimeY, Y12, label=stra, linewidth=3, linecolor=:black, linestyle=:solid)
        if L==3
            plot!(TimeY, Y13, label=stru, linewidth=2.5, linecolor=:black, linestyle=:dot)
        end
    end
    # Axis settings
    xlims!(FirstY, LastY)
    ylims!(0, maxD)
    xticks!(ceil(FirstY/Tgap)*Tgap:Tgap:floor(LastY/Tgap)*Tgap,[""])
    # Grid and title
    title!("Capacity investment (drilling), mmbbl/d per year")
    xgrid!(true)
    ygrid!(true)
    # Legend settings
    if L==1
        plot!(legend = false)
    else
        plot!(legend=:topright, legendfontsize=fsize)
    end
    # Font
    if serif==true
        plot!(tickfontsize=fsize, titlefontsize=fsize, fontfamily="Computer Modern")
    else
        plot!(tickfontsize=fsizess, titlefontsize=fsizess)
    end

    # Plot production
    qplot = plot(TimeY, Y21, label=strb, linewidth=1.5, linecolor=:red, linestyle=:solid)
    if L>=2
        plot!(TimeY, Y22, label=stra, linewidth=3, linecolor=:black, linestyle=:solid)
        if L==3
            plot!(TimeY, Y23, label=stru, linewidth=2.5, linecolor=:black, linestyle=:dot)
        end
    end
    # Axis settings
    xlims!(FirstY, LastY)
    ylims!(0, maxQ)
    xticks!(ceil(FirstY/Tgap)*Tgap:Tgap:floor(LastY/Tgap)*Tgap,[""])
    # Grid and title
    title!("Production, mmbbl/d")
    xgrid!(true)
    ygrid!(true)
    # Legend settings
    plot!(legend = false)
    # Font
    if serif==true
        plot!(tickfontsize=fsize, titlefontsize=fsize, fontfamily="Computer Modern")
    else
        plot!(tickfontsize=fsizess, titlefontsize=fsizess)
    end    

    # Plot oil price
    pplot = plot(TimeY, Y31, linewidth=1.5, linecolor=:red, linestyle=:solid)
    if L>=2
        plot!(TimeY, Y32, linewidth=3, linecolor=:black, linestyle=:solid)
        if L==3
            plot!(TimeY, Y33, linewidth=2.5, linecolor=:black, linestyle=:dot)
        end
    end
    # Axis settings
    xlims!(FirstY, LastY)
    ylims!(0, maxP)
    xticks!(ceil(FirstY/Tgap)*Tgap:Tgap:floor(LastY/Tgap)*Tgap)
    ytstring = ["\$0", "\$50", "\$100", "\$150", "\$200", "\$250"]
    yticks!(1:50:maxP,ytstring)
    # Grid and title
    title!("Oil price, \$/bbl")
    xgrid!(true)
    ygrid!(true)
    # Legend settings
    plot!(legend = false)
    # Font
    if serif==true
        plot!(tickfontsize=fsize, titlefontsize=fsize, fontfamily="Computer Modern")
    else
        plot!(tickfontsize=fsizess, titlefontsize=fsizess)
    end

    # Combine plots into a single figure
    if serif==true
        combined_plot = plot(dplot, qplot, pplot, 
            layout = @layout([a; b; c]), size=(1000, 1100))
    else
        combined_plot = plot(dplot, qplot, pplot, 
            layout = @layout([a; b; c]), size=(1000, 850))
    end
    return combined_plot
end



# ==============================================================
# Zoomed in plot for prices during the first few years
# ==============================================================
function priceplotzoomed(Oi::NamedTuple, Od::NamedTuple, Om::NamedTuple, 
    T::Int64, serif::Bool = false, L::Int64 = 3)
    # INPUTS
    # Oi: tuple of baseline outcomes
    # Od: tuple of anticipated demand decline outcomes
    # Om: tuple of myopic demand decline outcomes
    # T: number of periods per year
    # serif: true/false flag for whether to use a serif font
    # L: number of lines to plot (baseline, anticipated decline, unanticipated decline)

    # OUTPUTS
    # price_zoom_plot: plot of prices during first few Years

    # First set up time axis and region labels
    FirstY = 2023; LastY = 2035;        # first and last year to plot
    NTplot = (LastY-FirstY) * T + 1;  # number of times to plot
    TimeY = [FirstY:1/T:LastY];       # vector of time in calendar years
 
    # Extract price data
    Y31 = vec(sum(Oi.pvec[1:NTplot,:], dims=2));        # oil price
    Y32 = vec(sum(Od.pvec[1:NTplot,:], dims=2));  
    Y33 = vec(sum(Om.pvec[1:NTplot,:], dims=2)); 
    
    # More plot variables
    Tgap = 2;                           # gap between x-axis labels
    minP = 70; maxP = 90   # max Y axis values
    fsize = 15;                         # font size 
    fsizess = 15;                       # font size for slides

    # Legend labels
    strb = "baseline demand"
    stra = "anticipated demand decline"
    stru = "unanticipated demand decline"

    # Plot oil price
    price_zoom_plot = plot(TimeY, Y31, linewidth=1.5, linecolor=:red, linestyle=:solid)
    if L>=2
        plot!(TimeY, Y32, linewidth=3, linecolor=:black, linestyle=:solid)
        if L==3
            plot!(TimeY, Y33, linewidth=2.5, linecolor=:black, linestyle=:dot)
        end
    end
    # Axis settings
    xlims!(FirstY, LastY)
    ylims!(minP, maxP)
    xticks!(ceil(FirstY/Tgap)*Tgap:Tgap:floor(LastY/Tgap)*Tgap)
    ytstring = ["\$70", "\$75", "\$80", "\$85", "\$90"]
    yticks!(minP:5:maxP,ytstring)
    # Grid and title
    title!("Oil price, \$/bbl")
    xgrid!(true)
    ygrid!(true)
    # Legend settings
    plot!(legend = false)
    # Font
    if serif==true
        plot!(tickfontsize=fsize, titlefontsize=fsize, fontfamily="Computer Modern")
    else
        plot!(tickfontsize=fsizess, titlefontsize=fsizess)
    end

    return price_zoom_plot
end



# ==============================================================
# Function to plot production for a region
# ==============================================================
function regionplot(qveci::Matrix{Float64}, qvecd::Matrix{Float64}, qvecm::Matrix{Float64}, 
    T::Int64, i::Int64, serif::Bool = false, leg::Bool = false, yaxt::Bool = false, L::Int64 = 3)
    # INPUTS
    # qveci: matrix of production for baseline
    # qvecd: matrix of production for anticipated decline
    # qvecm: matrix of production for unanticipated decline
    # T: number of periods per year
    # i: region index
    # serif: true/false flag for whether to use a serif font
    # leg: true/false flag for whether to include legend
    # yaxt: true/false flag for whether to include y-axis title
    # L: number of lines to plot (baseline, anticipated decline, unanticipated decline)

    # OUTPUTS
    # rplot: plot of production for the region

    # First set up time axis and region labels
    FirstY = 2023; LastY = 2119;        # first and last year to plot
    NTplot = (LastY-FirstY) * T + 1;  # number of times to plot
    TimeY = [FirstY:1/T:LastY];       # vector of time in calendar years
    Tgap = 20;                          # gap between x-axis labels
    maxQ = 35;                          # max Y axis value
    fsize = 19;                         # font size 
    ftsize = 24;                        # title font size
    fsizess = 15;                       # font size for slides

    # Extract production for baseline, anticipated decline, and unanticipated decline
    Y1 = qveci[1:NTplot,i]      # baseline
    Y2 = qvecd[1:NTplot,i]      # anticipated decline
    Y3 = qvecm[1:NTplot,i]      # unanticipated decline 

    # Legend labels
    strb = "baseline demand"
    stra = "anticipated demand decline"
    stru = "unanticipated demand decline"

    # Plot production
    rplot = plot(TimeY, Y1, label=strb, linewidth=1.5, linecolor=:red, linestyle=:solid)
    if L>=2
        plot!(TimeY, Y2, label=stra, linewidth=3, linecolor=:black, linestyle=:solid)
        if L==3
            plot!(TimeY, Y3, label=stru, linewidth=2.5, linecolor=:black, linestyle=:dot)
        end
    end
    # Axis settings
    xlims!(FirstY, LastY)
    ylims!(0, maxQ)
    xticks!(ceil(FirstY/Tgap)*Tgap:Tgap:floor(LastY/Tgap)*Tgap)
    # Grid
    xgrid!(true)
    ygrid!(true)
    # y-axis title
    if yaxt==true
        ylabel!("Production, mmbbl/d")
    end
    # Legend settings
    if leg==true && L>=2
        if serif==true
            plot!(legend=:topright, legendfontsize=fsize)
        else
            plot!(legend=:topright, legendfontsize=fsizess)
        end
    else
        plot!(legend = false)
    end
    # Font
    if serif==true
        plot!(tickfontsize=fsize, titlefontsize=fsize, yguidefontsize=fsize, fontfamily="Computer Modern")
    else
        plot!(tickfontsize=fsize, titlefontsize=ftsize, yguidefontsize=fsize)
    end
    return rplot
end
