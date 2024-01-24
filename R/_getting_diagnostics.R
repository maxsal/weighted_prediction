which(phers_calculated[, name == "t0_pwide_univariable_results"])

ukb_phers_list

tmp <- ukb_phers_list[[which(phers_calculated[, name == "t0_pwide_univariable_results"])]]

#####
.clean_auc <- function(data, outcome = "case", exposure = "phers") {
    tmp <- suppressMessages(pROC::roc(
        response  = data[[outcome]],
        predictor = data[[exposure]],
        family    = binomial(),
        ci        = TRUE
    ))

    output <- data.frame(
        auc = tmp$auc,
        lower = tmp$ci[1],
        upper = tmp$ci[3]
    )
    output$auc_print <- paste0(
        round(output$auc, 3),
        " (",
        round(output$lower, 3),
        ", ",
        round(output$upper, 3),
        ")"
    )
    return(output)
}

getTopEffects <- function(prob, data, outcome = "case", exposure, covs = NULL, pl = FALSE) {
    riskBin              <- paste0("Top_", prob)
    data[[riskBin]] <- ifelse(
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = 1 - prob),
        1, 0)
    if (!is.null(covs)) {
        tmp_f <- paste0(outcome, " ~ ", riskBin, " + ", paste(covs, collapse = " + "))
    } else {
        tmp_f <- paste0(outcome, " ~ ", riskBin)
    }
    vars <- c(outcome, riskBin)
    tmp_data <- na.omit(data[, ..vars])
    if (!is.null(covs)) {
        vars <- c(outcome, riskBin, covs)
    } else {
        vars <- c(outcome, riskBin)
    }
    fitAsso2 <- logistf::logistf(tmp_f, data = tmp_data, plcontrol = logistpl.control(maxit = 1E5), pl = pl)
    getValues(fitAsso2, riskBin)
}

getTopPower <- function(prob, data, outcome, exposure) {
    riskBin <- glue("Top_{prob}")
    data[[riskBin]] <- ifelse(
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = 1 - prob),
        1, 0)
    ptable <- table("Top" = data[[riskBin]], "trait" = data[[outcome]])
    if (nrow(ptable) != 2) {
        return(c("power" = NA, "h" = NA))
    }
    ntop <- sum(ptable[, 2])
    nrest <- sum(ptable[, 1])
    ptop <- ptable[2, 2] / ntop
    prest <- ptable[2, 1] / nrest
    h <- pwr::ES.h(ptop, prest)
    unlist(pwr::pwr.2p2n.test(h = h, n1 = ntop, n2 = nrest, sig.level = 0.05, power = NULL)[c("power", "h")])
}

log10toP <- function(log10P) {
    log10P <- abs(as.numeric(log10P))
    if (is.na(log10P)) {
        return(NA)
    }
    if (log10P > 300) {
        part1 <- log10P %/% 100 * 100
        part2 <- log10P - part1
        if (part2 != 0) {
            P <- format(signif(10^-part2, 3), scientific = T)
            P <- paste(as.numeric(gsub("e-.+", "", P)), "e-", as.numeric(gsub(".+-", "", P), sep = "") + part1, sep = "")
        } else {
            P <- paste("1e-", part1, sep = "")
        }
    } else {
        P <- signif(10^-log10P, 3)
    }
    return(as.character(P))
}

getValues <- function(fitAsso, exposure) {
    sebeta <- sqrt(diag(vcov(fitAsso)))
    names(sebeta) <- fitAsso$terms
    BETA <- round(fitAsso$coefficient[exposure], 4)
    SEBETA <- round(sebeta[exposure], 4)

    vars <- diag(fitAsso$var)
    names(vars) <- fitAsso$terms
    LOG10P <- -pchisq((fitAsso$coefficient[exposure]^2 / vars[exposure]), 1, lower.tail = F, log = T) / log(10)
    P <- log10toP(LOG10P)

    OR <- round(exp(fitAsso$coefficient[exposure]), 4)
    CI1 <- round(exp(fitAsso$ci.lower[exposure]), 4)
    CI2 <- round(exp(fitAsso$ci.upper[exposure]), 4)
    out <- c(BETA, SEBETA, P, OR, CI1, CI2, signif(LOG10P, 4),
        paste0(
            format(round(OR, 2), nsmall = 2, big.mark = ","),
            " (",
            format(round(CI1, 2), nsmall = 2, big.mark = ","),
            ", ",
            format(round(CI2, 2), nsmall = 2, big.mark = ","),
            ")")
        )
    names(out) <- paste0(exposure, "_", c("beta", "se_beta", "p_value", "or", "lower", "upper", "log10p", "print"))
    return(out)
}
#####

get_diagnostics <- function(data, outcome, exposure, covs = NULL, pctile_or = TRUE, pl = FALSE) {
    out <- list()
    # formulas
        # without covs
        f_nocovs <- as.formula(paste0(outcome, " ~ ", exposure))
        # with covs
        if (!is.null(covs)) {
            f_covs <- paste0(outcome, " ~ ", exposure, " + ", paste(covs, collapse = " + "))
            f_covs <- as.formula(f_covs)
        } else {
            f_covs <- f_nocovs
        }

    # OR (glm)
    glm_mod <- glm(f_covs, data = data, family = binomial())
    out$beta    <- glm_mod$coefficient[exposure]
    out$se_beta <- sqrt(diag(vcov(glm_mod))[exposure])
    out$or      <- exp(out$beta)
    glm_mod_confint <- suppressMessages(confint(glm_mod))
    out$or_lower <- exp(glm_mod_confint[exposure, 1])
    out$or_upper <- exp(glm_mod_confint[exposure, 2])
    out$or_print <- paste0(
        format(round(out$or, 2), nsmall = 2),
        " (",
        format(round(out$or_lower, 2), nsmall = 2),
        ", ",
        format(round(out$or_upper, 2), nsmall = 2),
        ")"
    )
    out$log10p <- -log10(summary(glm_mod)$coefficients[exposure, 4])

    # AUC (pROC::roc)
    tmp_auc <- suppressMessages(pROC::roc(
            response  = data[[outcome]],
            predictor = data[[exposure]],
            family    = binomial(),
            ci        = TRUE
    ))
    out$auc       <- tmp_auc$auc[1]
    out$auc_lower <- tmp_auc$ci[1]
    out$auc_upper <- tmp_auc$ci[3]
    out$auc_print <- paste0(
        format(round(out$auc, 3), nsmall = 3),
        " (",
        format(round(out$auc_lower, 3), nsmall = 3),
        ", ",
        format(round(out$auc_upper, 3), nsmall = 3),
        ")"
    )

    glm_no_cov_mod <- glm(f_nocovs, data = as.data.frame(data), family = binomial())
    # R2 (rcompanion::nagelkerke)
    nagel_out <- rcompanion::nagelkerke(fit = glm_no_cov_mod)
    out[["R2 (McFadden)"]]                  <- nagel_out$Pseudo.R.squared.for.model.vs.null[1]
    out$"R2 (Cox and Snell [ML])"           <- nagel_out$Pseudo.R.squared.for.model.vs.null[2]
    out$"R2 (Nagelkerke [Cragg and Uhler])" <- nagel_out$Pseudo.R.squared.for.model.vs.null[3]

    # HL test (ResourceSelection::hoslem.test)
    hl_out <- ResourceSelection::hoslem.test(glm_no_cov_mod$y, fitted(glm_no_cov_mod), g = 10)
    out$HosmerLemeshow_ChiSq <- hl_out$statistic[1]
    out$HosmerLemeshow_P     <- hl_out$p.value
    out$hl_print <- paste0(
        format(round(out$HosmerLemeshow_ChiSq, 2), nsmall = 2),
        " (",
        formatC(out$HosmerLemeshow_P, format = "e", digits = 2),
        ")"
    )
    # Brier score (DescTools)
    out$Brier_score <- DescTools::BrierScore(glm_no_cov_mod)
    out$bs_print <- format(round(out$Brier_score, 5), nsmall = 5)

    # %ile-based OR
    if (pctile_or) {
        probs <- c(0.01, 0.02, 0.05, 0.1, 0.25)
        for (pr in probs) {
            out <- c(out, getTopEffects(prob = pr, data = data, outcome = outcome, exposure = exposure, covs = covs, pl = pl))
        }
        
        # determine percentile with at least 80% power (or the most powered percentiles)
        ptest <- sapply(seq(0.005, 0.5, by = 0.005), \(x) getTopPower(x, data = data, outcome = outcome, exposure = exposure))
        presults <- data.table("Top" = seq(0.005, 0.5, by = 0.005), t(ptest))
        #  presults <- data.table('Top' = seq(0.005, 0.5, by=0.005), t(ptest))
        maxp <- presults$Top[which(presults$h == max(presults$h, na.rm = TRUE))[1]]

        p80 <- presults$Top[which(presults$power >= 0.80)]
        underpowered <- F
        if (length(p80) == 0) {
            p80 <- maxp
            underpowered <- T
        } else {
            p80 <- p80[1]
        }

        prs.p80 <- getTopEffects(prob = p80, data = data, outcome = outcome, exposure = exposure, covs = covs, pl = pl)
        names(prs.p80) <- gsub(p80, "MinPower80", names(prs.p80))

        out <- c(out, "MinPower80" = p80, prs.p80, "Top_Underpowered" = underpowered)
    } 

    data.table(
        stat  = names(out),
        value = unlist(out)
    )
}

aimTwo::get_bin_diagnostics(data = tmp, outcome = "case", exposure = "phers")



