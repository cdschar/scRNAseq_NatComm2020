library("beeswarm")

gType = "both" #bee,bar,both

gHeight = 2.5;
gWidth = 2.5
cex.main = 1.25;
cex.lab = 0.75
cex.axis = 0.75;
font.lab = 2;
font.axis = 1;
line.width = 0.2
mean.width = 0.2;
sd.width = 0.1;
gMgp = c(1.2, 0.5, 0)
gMai = c(0.7, 0.4, 0.2, 0.1);
gBty = "o"
mean.lwd = 2
sd.lwd = 1
gPch = 19
gLas = 3
gBuffer = 0.1

bbPlot = function(data, groups, type = c(gType), outFile = c(NA), cols, height = c(gHeight), width = c(gWidth), main = c(NA), xlab = c(NA), ylab = c(NA), ylim = c(NA), cex.main = c(1.25), cex.lab = c(0.75), cex.axis = c(0.75), las= c(gLas), line.lwd = c(1), line.width = c(0.2), mai = c(gMai), mgp = c(gMgp), bty = c(gBty), font.lab = c(2), font.axis = c(2), cex.dots = c(1), pch = c(gPch), yMaxBuffer = c(gBuffer), names = c(NA), ...) {

	#output to pdf or to X11
	if (!is.na(outFile)) cairo_pdf(outFile, height = height, width = width) else X11(height = height, width = width)
	
	#set par options
	par(mai = mai, font = font.lab, mgp = mgp, las = las, bty = bty, family = "Arial");

	#set colors if not provided
	if (all(is.na(cols))) cols = rgb(0.5, 0.5, 0.5)

	#calculate group average and variation
	means = aggregate(as.numeric(data) ~ groups, FUN = "mean", na.rm = T)
	var = aggregate(as.numeric(data) ~ groups, FUN = "var", na.rm = T)
	sd = var; sd[, 2] = sqrt(sd[, 2])
	
	#set sample names
	if(all(is.na(names))) names = unique(groups);

	#output beeswarm over boxplot
	if(type == "both") {
	
		if (all(is.na(ylim))) {
		
			yMax = max(as.numeric(data)) + yMaxBuffer * max(as.numeric(data))
			#yMax = (max(means[, 2], na.rm = TRUE) + max(sd[, 2], na.rm = TRUE)) * buffer;
			
			ylim = c(0, yMax)
			#ylim = c(0, ceiling(yMax / (10 ^ (nchar(round(yMax)) - 1))) * 10 ^ (nchar(round(yMax)) - 1));
		}

		#output bar plot
		bp = barplot(means[, 2], names = names, main = main, cex.main = cex.main, cex.lab = cex.lab, 
					 cex.axis = cex.axis, cex.name = cex.lab, cex.lab = cex.lab, ylim = ylim, col = "white", border = cols, 
					 ylab = ylab, font.lab = font.lab, font.axis = font.axis, mgp = mgp, las = las, ...);
	
		#overlay beeswarm using bp[,1] as the x axis coordinates
		beeswarm(data ~ groups, ylim=ylim, add=T, col=cols, at=bp[,1], cex = cex.dots, pch = pch, ...)
	
		for (i in 1:dim(means)[1]) {

			lines(c(bp[i, 1], bp[i, 1]), c(means[i, 2], means[i, 2] + sd[i, 2]), lwd = line.lwd, col = cols[i]);
			lines(c(bp[i, 1] - line.width, bp[i, 1] + line.width), c(means[i, 2] + sd[i, 2], means[i, 2] + sd[i, 2]), lwd = line.lwd, col = cols[i]);
		}
		
		if (!is.na(outFile)) dev.off();
	}

	#output beeswarm only
	if(type == "bee") {
	
		ymin = min(as.numeric(data)) - 0.1 * min(as.numeric(data))
		ymax = max(as.numeric(data)) + 0.1 * max(as.numeric(data))
		if (all(is.na(ylim))) ylim = c(ymin, ymax)

		b = beeswarm(data ~ groups, labels = names, cex.lab = cex.lab, cex.axis = cex.axis, col = cols, ylab = ylab, xlab = xlab, 
			main = main, ylim = ylim, cex = cex.dots, pch = pch, las = las, bty = bty, ...) #

		for (i in 1:length(unique(groups))) {

			lines(c(i - mean.width, i + mean.width), rep(means[i, 2], 2), lwd = mean.lwd)
			lines(rep(i, 2), c(means[i, 2] - sd[i, 2], means[i, 2] + sd[i, 2]), lwd = sd.lwd)
			lines(c(i - sd.width, i + sd.width), c(means[i, 2] + sd[i, 2], means[i, 2] + sd[i, 2]), lwd = sd.lwd)
			lines(c(i - sd.width, i + sd.width), c(means[i, 2] - sd[i, 2], means[i, 2] - sd[i, 2]), lwd = sd.lwd)
		}

		if (!is.na(outFile)) dev.off();
	
	}

	#output barplot only
	if(type == "bar") {

		if (all(is.na(ylim))) {
		
			yMax = max(as.numeric(data)) + yMaxBuffer * max(as.numeric(data))
			#yMax = (max(means[, 2], na.rm = TRUE) + max(sd[, 2], na.rm = TRUE)) * buffer;
			
			ylim = c(0, yMax)
			#ylim = c(0, ceiling(yMax / (10 ^ (nchar(round(yMax)) - 1))) * 10 ^ (nchar(round(yMax)) - 1));
		}

		bp = barplot(means[, 2], names = names, main = main, cex.main = cex.main, cex.lab = cex.lab, 
					 cex.axis = cex.axis, cex.name = cex.lab, cex.lab = cex.lab, ylim = ylim, col = cols, border = cols, 
					 ylab = ylab, font.lab = font.lab, font.axis = font.axis, mgp = mgp, las = las, ...);
					 
		for (i in 1:dim(means)[1]) {

			lines(c(bp[i, 1], bp[i, 1]), c(means[i, 2], means[i, 2] + sd[i, 2]), lwd = line.lwd, col = cols[i]);
			lines(c(bp[i, 1] - line.width, bp[i, 1] + line.width), c(means[i, 2] + sd[i, 2], means[i, 2] + sd[i, 2]), lwd = line.lwd, col = cols[i]);
		}

		if (!is.na(outFile)) dev.off();
	}

}
