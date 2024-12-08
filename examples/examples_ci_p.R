# Examples for single values
ci_p(x=5, n=15, intervalType="wald", conf.level=0.90)
ci_p(x=7, n=17, intervalType="wald", conf.level=0.95)
ci_p(x=9, n=19, intervalType="wald", conf.level=0.97)

# Examples for multiple values
ci_p(x=c(5, 7, 9),
     n=c(15, 17, 19),
     intervalType="wald",
     conf.level=c(0.90, 0.95, 0.97))
