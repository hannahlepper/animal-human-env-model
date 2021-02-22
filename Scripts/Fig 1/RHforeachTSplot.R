#sort out julia's list
R <- unlist(R)
#For testing - dummy data
#R = rbeta(3 * 5, 1, 1)
Rdf <- data.frame(
        R = R,
        R_type = factor(rep(c("R[H]", "R[A]", "R[E]"), times = 5)),
        TS = rep(c("Baseline", "Balanced", "Environment\ndominated", "Animal\ndominated", "Human\ndominated"), each = 3))

png("plots/barplotFig1B.png", width = 8, height = 7, units = "cm", res = 300)
p <- ggplot(Rdf, aes(TS, R, fill = R_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_discrete("", labels = parse(text = levels(Rdf$R_type))) +
    labs(x = "", y = "Fraction of population carrying\n resistant bacteria") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
print(p)
dev.off()
