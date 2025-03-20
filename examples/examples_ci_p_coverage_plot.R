ci_p_coverage_plot(n=10,
                   intervalType="wald",
                   conf.level=0.95,
                   ylim=c(0.8, 1), las=1)

ci_p_coverage_plot(n=10,
                   intervalType="clopper_pearson",
                   conf.level=0.95,
                   ylim=c(0.9, 1), las=1)

ci_p_coverage_plot(n=10,
                   intervalType="wilson",
                   conf.level=0.95,
                   ylim=c(0.9, 1), las=1)

ci_p_coverage_plot(n=10,
                   intervalType="jeffreys",
                   conf.level=0.95,
                   ylim=c(0.9, 1), las=1)

ci_p_coverage_plot(n=10,
                   intervalType="agresti_coull",
                   conf.level=0.95,
                   ylim=c(0.9, 1), las=1)

ci_p_coverage_plot(n=10,
                   intervalType="arcsine",
                   conf.level=0.95,
                   ylim=c(0.9, 1), las=1)
