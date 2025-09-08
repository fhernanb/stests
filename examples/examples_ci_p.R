# Examples for alternative IC
ci_p(x=5, n=15, intervalType="wald", conf.level=0.90)
ci_p(x=5, n=15, intervalType="agresti_coull", conf.level=0.90)
ci_p(x=5, n=15, intervalType="rindskopf", conf.level=0.90)
ci_p(x=5, n=15, intervalType="clopper_pearson", conf.level=0.90)
ci_p(x=5, n=15, intervalType="arcsine_anscombe", conf.level=0.90)
ci_p(x=5, n=15, intervalType="arcsine", conf.level=0.90)

# Examples for multiple values
ci_p(x=c(5, 7, 9),
     n=c(15, 17, 19),
     intervalType="wald",
     conf.level=c(0.90, 0.95, 0.97))
