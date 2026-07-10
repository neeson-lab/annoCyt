.msg = function(text, verbose = FALSE, level = 0, value = NULL) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  indent = strrep("  ", level)
  if (is.null(value)) message(indent, text)
  else                message(indent, text, ": ", value)
}


# ---- Internal: extract marker name from a rule string ----------------------
.extractMarker = function(rule, prefix, direction, fallback_marker) {
  tag = paste0(prefix, "\\.", direction)
  if (rule == paste0(prefix, ".", direction)) {
    return(fallback_marker)
  }
  sub(paste0("(.*)", tag, "$"), "\\1", rule)
}


# ---- Internal: vectorised threshold comparison ------------------------------
.evalThresholdRules = function(rules_df, prefix, fallback_marker, get_val_fn) {

  rules      = rules_df$rule
  thresholds = rules_df$threshold

  direction = ifelse(grepl(paste0(prefix, "\\.greaterthan"), rules),
                      "greaterthan", "lessthan")

  markers = mapply(.extractMarker,
                    rule            = rules,
                    direction       = direction,
                    MoreArgs        = list(prefix          = prefix,
                                           fallback_marker  = fallback_marker),
                    SIMPLIFY        = TRUE)

  vals = vapply(markers, get_val_fn, numeric(1))

  pass = ifelse(
    direction == "greaterthan",
    ifelse(thresholds <= 0, vals > thresholds, vals >= thresholds),
    ifelse(thresholds <= 0, vals <= thresholds, vals < thresholds)
  )

  if (all(pass)) "pass" else "fail"
}


# ---- .ruleExpr --------------------------------------------------------------
.ruleExpr = function(filtered_lt, expression_data, filtered_rt) {

  rt = filtered_rt[grepl("e\\.", filtered_rt$rule), , drop = FALSE]
  if (nrow(rt) == 0L) return("pass")

  raw_marker = filtered_lt$value[[1L]][[1L]]

  get_expr = function(marker) {
    col_idx = match(marker, colnames(expression_data))
    if (is.na(col_idx)) stop("Marker not found in expression data: ", marker)
    as.double(expression_data[[col_idx]][1L])
  }

  gt_rt = rt[grepl("e\\.greaterthan", rt$rule), , drop = FALSE]
  lt_rt = rt[grepl("e\\.lessthan",    rt$rule), , drop = FALSE]

  results = character(0L)

  if (nrow(gt_rt) > 0L) {
    results = c(results, .evalThresholdRules(gt_rt, "e", raw_marker, get_expr))
  }

  if (nrow(lt_rt) > 0L) {
    lt_fallback = gsub("neg", "", raw_marker)
    results = c(results, .evalThresholdRules(lt_rt, "e", lt_fallback, get_expr))
  }

  if (any(results == "fail")) "fail" else "pass"
}


# ---- .ruleScore -------------------------------------------------------------
.ruleScore = function(filtered_lt, ranking_scores, filtered_rt) {

  rt = filtered_rt[grepl("s\\.", filtered_rt$rule), , drop = FALSE]
  if (nrow(rt) == 0L) return("pass")

  default_label = filtered_lt$label[[1L]]

  get_score = function(marker) {
    col_name = paste0(marker, "_TotalScore")
    col_idx  = match(col_name, names(ranking_scores))
    if (is.na(col_idx)) stop("Score column not found: ", col_name)
    as.double(ranking_scores[[col_idx]][1L])
  }

  .evalThresholdRules(rt, "s", default_label, get_score)
}


# ---- .ruleMax ---------------------------------------------------------------
.ruleMax = function(rule_x, possible_label, ranking_scores, verbose) {

  children_str = regmatches(rule_x,
                             regexpr("(?<=\\()[^()]+(?=\\))", rule_x, perl = TRUE))

  children = trimws(strsplit(children_str, ",")[[1L]])
  col_idx = match(paste0(children, "_TotalScore"), names(ranking_scores))
  scores  = unlist(ranking_scores[1L, col_idx])
  names(scores) = children
  winners = names(scores)[scores == max(scores)]
  if (length(winners) == 1L && possible_label == winners) {
    .msg("Passed: max rule", verbose, level = 4L)
    return(winners)
  } else if (length(winners) == 1L && possible_label != winners) {
    .msg("Failed: max rule", verbose, level = 4L)
    return(winners)
  } else if (length(winners) > 1L && possible_label %in% winners){
    msgw = paste(winners, sep = ", ")
    .msg("Warning: multiple cell types have same score. Both may be assigned unless other rules are in place. Cell types with same score: ",
         verbose, level = 4L, value = msgw)
    return(possible_label)
  } else if (!(possible_label %in% winners)) {
    .msg("Failed: max rule", verbose, level = 4L)
    return("empty")
  }
  }




# ---- Shared rule-evaluation core --------------------------------------------
# full_rulesdf: complete rules table, needed to look up rules for sibling
#               labels referenced by not. and isct. rules.
.evalLabelRules = function(lineagedf,
                           expression_data,
                           ranking_scores,
                           cluster_label_full,
                           rulesdf,
                           full_rulesdf,
                           verbose) {
  check = "nil"

  # The shared parent of siblings at this level, derived from the lineage rows
  # for the current label (already filtered to one label before calling).
  parent = lineagedf$parent[[1L]]

  for (row_x in seq_len(nrow(rulesdf))) {
    rule_x         = rulesdf$rule[[row_x]]
    possible_label = rulesdf$label[[1L]]

    if (grepl("max", rule_x, fixed = TRUE)) {
      label_max = .ruleMax(rule_x, possible_label, ranking_scores, verbose)
      check = if (label_max != possible_label) "fail"
      else if (check != "fail") "pass" else check

    } else if (grepl("s\\.", rule_x)) {
      ans = .ruleScore(
        lineagedf[lineagedf$label == possible_label, , drop = FALSE],
        ranking_scores,
        rulesdf[row_x, , drop = FALSE]
      )
      check = if (ans == "fail") "fail" else if (check != "fail") "pass" else check

    } else if (grepl("e\\.", rule_x)) {
      ans = .ruleExpr(
        lineagedf[lineagedf$label == possible_label, , drop = FALSE],
        expression_data,
        rulesdf[row_x, , drop = FALSE]
      )
      check = if (ans == "fail") "fail" else if (check != "fail") "pass" else check

    } else if (grepl("not\\.", rule_x, fixed = TRUE)) {
      ct_x   = sub("not\\.", "", rule_x)
      # Sibling: same parent as current label, different label name
      sub_lin   = lineagedf[lineagedf$label == ct_x &
                               lineagedf$parent == sibling_parent, , drop = FALSE]
      sub_rules = full_rulesdf[full_rulesdf$label == ct_x, , drop = FALSE]
      ifct = .ruleIsLabel(sub_lin, expression_data, ranking_scores,
                            cluster_label_full, sub_rules, full_rulesdf)
      check = if (ifct == "not.celltype" && check != "fail") "pass"
      else if (ifct == "is.celltype") "fail" else check

    } else if (grepl("isct\\.", rule_x, fixed = TRUE)) {
      ct_x   = sub("isct\\.", "", rule_x)
      # Sibling: same parent as current label, different label name
      sub_lin   = lineagedf[lineagedf$label == ct_x &
                               lineagedf$parent == sibling_parent, , drop = FALSE]
      sub_rules = full_rulesdf[full_rulesdf$label == ct_x, , drop = FALSE]
      ifct = .ruleIsLabel(sub_lin, expression_data, ranking_scores,
                           cluster_label_full, sub_rules, full_rulesdf)
      check = if (ifct == "is.celltype" && check != "fail") "pass"
      else if (ifct == "not.celltype") "fail" else check
    }
  }

  if (check == "pass") "is.celltype" else "not.celltype"
}


# ---- .ruleIsLabel / .ruleNotLabel -------------------------------------------
.ruleIsLabel = function(lineagedf, expression_data, ranking_scores,
                        cluster_label_full, rulesdf, full_rulesdf) {
  .evalLabelRules(lineagedf, expression_data, ranking_scores,
                  cluster_label_full, rulesdf, full_rulesdf, verbose = FALSE)
}

# .ruleNotLabel = function(lineagedf, expression_data, ranking_scores,
#                          cluster_label_full, rulesdf, full_rulesdf) {
#   .evalLabelRules(lineagedf, expression_data, ranking_scores,
#                   cluster_label_full, rulesdf, full_rulesdf, verbose = FALSE)
# }


# ---- .runRules --------------------------------------------------------------
.runRules = function(lineagedf,
                     expression_data,
                     ranking_scores,
                     cluster_label_full,
                     rulesdf,
                     verbose) {

  parent_label = if (length(cluster_label_full) == 0L) "root"
  else cluster_label_full[length(cluster_label_full)]

  filtered_lt = lineagedf[lineagedf$parent == parent_label, , drop = FALSE]

  if (nrow(filtered_lt) == 0L) {
    # message("Error: no children found for parent '", parent_label, "'")
    return(c(cluster_label_full, "stop"))
  }

  possible_labels = unique(filtered_lt$label)

  # Named status vector — stores "pass"/"fail"/"nil" per candidate label
  status = stats::setNames(rep("nil", length(possible_labels)), possible_labels)

  for (i in possible_labels) {

    if (status[[i]] != "nil") next # already been evaluated in not/is rule, skip
    .msg("Checking rules for possible label", verbose, level = 3L, value = i)

    rt_i  = rulesdf[rulesdf$label == i, , drop = FALSE]
    lt_i  = filtered_lt[filtered_lt$label == i, , drop = FALSE]
    check = "nil"

    for (row_x in seq_len(nrow(rt_i))) {
      rule_x = rt_i$rule[[row_x]]

      if (grepl("max", rule_x, fixed = TRUE)) {
        label_max = .ruleMax(rule_x, i, ranking_scores, verbose)
        check = if (label_max != i) "fail"
        else if (check != "fail") "pass" else check

      } else if (grepl("s\\.", rule_x)) {
        ans = .ruleScore(lt_i, ranking_scores, rt_i[row_x, , drop = FALSE])
        if (ans == "fail") .msg("Failed: compare singscore to threshold", verbose, level = 4L, value = rule_x)
        else               .msg("Passed: compare singscore to threshold", verbose, level = 4L, value = rule_x)
        check = if (ans == "fail") "fail" else if (check != "fail") "pass" else check

      } else if (grepl("e\\.", rule_x)) {
        ans = .ruleExpr(lt_i, expression_data, rt_i[row_x, , drop = FALSE])
        if (ans == "fail") .msg("Failed: compare expression to threshold", verbose, level = 4L)
        check = if (ans == "fail") "fail" else if (check != "fail") "pass" else check

      } else if (grepl("not\\.", rule_x)) {
        ct_x = sub("not\\.", "", rule_x)

        if (ct_x %in% names(status)) {
          # Cache hit: status stores "pass"/"fail" for whether ct_x IS a cell
          # type. not. inverts: if ct_x passed -> current fails; if failed -> passes.
          cached = status[[ct_x]]

          if (cached == "nil"){
            eval = .ruleIsLabel(lt_i, expression_data, ranking_scores, cluster_label_full, rt_i, rulesdf)
            status[[ct_x]] = if (eval == "is.celltype") "pass" else if (eval == "not.celltype") "fail"
            cached = status[[ct_x]]
          }

          check  = if (cached == "pass") "fail"
          else if (cached == "fail" && check != "fail") "pass" else check
         } # else {
        #   # Fresh evaluation: ct_x is a sibling (same parent_label)
        #   sub_lin   = lineagedf[lineagedf$label == ct_x &
        #                            lineagedf$parent == parent_label, , drop = FALSE]
        #   sub_rules = rulesdf[rulesdf$label == ct_x, , drop = FALSE]
        #   ifct = .ruleNotLabel(sub_lin, expression_data, ranking_scores,
        #                         cluster_label_full, sub_rules, rulesdf)
        #   if (ifct == "not.celltype") .msg("Passed: not", verbose, level = 4L, value = ct_x)
        #   else                        .msg("Failed: not", verbose, level = 4L, value = ct_x)
        #   check = if (ifct == "not.celltype" && check != "fail") "pass"
        #   else if (ifct == "is.celltype") "fail" else check
        # }

      } else if (grepl("isct\\.", rule_x)) {
        ct_x = sub("isct\\.", "", rule_x)


        if (ct_x %in% names(status)) {
          # Cache hit: isct. needs ct_x to have passed.
          cached = status[[ct_x]]
          if (cached == "nil"){
            eval = .ruleIsLabel(lt_i, expression_data, ranking_scores, cluster_label_full, rt_i, rulesdf)
            status[[ct_x]] = if (eval == "is.celltype") "pass" else if (eval == "not.celltype") "fail"
            cached = status[[ct_x]]
          }

          check  = if (cached == "pass" && check != "fail") "pass"
          else if (cached == "fail") "fail" else check
         } # else {
        #   # Fresh evaluation: ct_x is a sibling (same parent_label)
        #   sub_lin   = lineagedf[lineagedf$label == ct_x &
        #                            lineagedf$parent == parent_label, , drop = FALSE]
        #   sub_rules = rulesdf[rulesdf$label == ct_x, , drop = FALSE]
        #   ifct = .ruleIsLabel(sub_lin, expression_data, ranking_scores,
        #                        cluster_label_full, sub_rules, rulesdf)
        #   if (ifct == "is.celltype") .msg("Passed: isct", verbose, level = 4L, value = ct_x)
        #   else                       .msg("Failed: isct", verbose, level = 4L, value = ct_x)
        #   check = if (ifct == "is.celltype" && check != "fail") "pass"
        #   else if (ifct == "not.celltype") "fail" else check
        # }
      }
    }

    status[[i]] = check
  }

  # isct. rule: when both the sibling label and sub-label pass, drop the sibling
  isct_rules = rulesdf[grepl("isct\\.", rulesdf$rule), , drop = FALSE]
  if (nrow(isct_rules) > 0L) {
    for (k in seq_len(nrow(isct_rules))) {
      main_lbl = sub("isct\\.", "", isct_rules$rule[[k]])
      sub_lbl  = isct_rules$label[[k]]
      if (main_lbl %in% names(status) &&
          status[[main_lbl]] == "pass" &&
          sub_lbl  %in% names(status) &&
          status[[sub_lbl]]  == "pass") {
        status = status[names(status) != main_lbl]
      }
   }
  }

  passed = names(status)[status == "pass"]
  n_pass = length(passed)

  if (n_pass == 1L) {
    .msg("All rules passed for", verbose, level = 3L, value = passed)
    return(c(cluster_label_full, passed))
  } else if (n_pass == 0L) {
    .msg("No child label passed all rules; unassigned", verbose, level = 3L)
    return(c(cluster_label_full, "stop"))
  } else {
    .msg("Warning: multiple children passed all their rules; both labels assigned.",
         verbose, level = 3L)
    return(c(cluster_label_full, paste(passed, collapse = ".")))
  }
}


# ---- .setLabelLevel ---------------------------------------------------------
.setLabelLevel = function(lineagedf,
                          expression_data,
                          ranking_scores,
                          cluster_label_full,
                          rulesdf,
                          verbose) {

  repeat {
    last = if (length(cluster_label_full) == 0L) "root"
    else cluster_label_full[length(cluster_label_full)]

    if (last == "stop" || !last %in% lineagedf$parent) break

    cluster_label_full = .runRules(lineagedf, expression_data, ranking_scores,
                                    cluster_label_full, rulesdf, verbose)
  }

  .msg("Finished annotating cluster. Labelled as", verbose, level = 1L,
       value = paste(cluster_label_full, collapse = "-"))
  cluster_label_full
}


#' Run annotation on SCE object
#' 
#' Determines the cell type label for every cluster using the lineage table, rules table, marker expression per cluster and singscore values. 
#' If a cluster passes one level of the hierarchy, it will attempt to annotate at the next level in the hierarchy until no more child cell types exist, or no child cell subtypes pass any rules. 
#' In the latter case, the cell type will be labelled as "unknown" at that particular level.
#' Can be subsetted to annotate only specified clusters. 
#' 
#' @details Set verbose = TRUE to determine which rules or clusters are not behaving as expected.
#'
#' @param expression_data  DataFrame of expression per marker per cluster; first column must be `cluster_id`,
#'                         remaining columns are marker expression values.
#' @param ranking_scores   DataFrame of singscore output per cluster; first column must be `cluster_id`,
#'                         remaining columns are `<label>_TotalScore` columns.
#' @param rulesdf          DataFrame of rules specifying how to call cell types. Must contain columns `label`, `rule`, `threshold`.
#' @param lineagedf        DataFrame of cell type hierarchy. Must contain columns `label`, `parent`, `value`.
#' @param which_clusters   character vector of cluster IDs to annotate,
#'                         or `"all"` (default).
#' @param verbose          logical; print progress messages (default `FALSE`).
#' @return A DataFrame of hierarchical annotation label per cluster
#' @export
#' 
#' @examples
#' \dontrun{
#' annoTbl = annotateClusters(eg_expr, eg_singscores, eg_rulesTable, eg_lineageTable)}
#' 
annotateClusters = function(expression_data,
                            ranking_scores,
                            rulesdf,
                            lineagedf,
                            which_clusters = "all",
                            verbose = FALSE) {

  suppressMessages({
    checkLineageTable(lineagedf, expression_data)
    checkRulesTable(rulesdf = rulesdf, lineagedf = lineagedf, exprdf = expression_data)
  })

  lineagedf = .updateLineageTable(lineagedf)

  if (!identical(which_clusters, "all")) {
    expression_data = expression_data[expression_data$cluster_id %in% which_clusters, , drop = FALSE]
    ranking_scores  = ranking_scores[ranking_scores$cluster_id   %in% which_clusters, , drop = FALSE]
  }

  expr_split = split(expression_data, expression_data$cluster_id)
  rank_split = split(ranking_scores,  ranking_scores$cluster_id)

  cluster_results = lapply(names(expr_split), function(j) {
    .msg("Annotating cluster", verbose, level = 0L, value = j)

    values = .setLabelLevel(
      lineagedf,
      expr_split[[j]],
      rank_split[[j]],
      character(0L),
      rulesdf,
      verbose
    )

    label_df = data.frame(
      cluster_id = j,
      t(stats::setNames(values, paste0("label", seq_along(values)))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    label_df
  })

  cluster_assignment = dplyr::bind_rows(cluster_results)

  # ---- Post-process labels --------------------------------------------------
  label_cols  = grep("^label", names(cluster_assignment), value = TRUE)
  parents_set = lineagedf$parent

  for (col in label_cols[-1L]) {
    prev_col       = label_cols[match(col, label_cols) - 1L]
    is_stop        = !is.na(cluster_assignment[[col]]) & cluster_assignment[[col]] == "stop"
    prev_is_parent = cluster_assignment[[prev_col]] %in% parents_set
    cluster_assignment[[col]] = ifelse(is_stop & prev_is_parent,
                                        "unknown",
                                        cluster_assignment[[col]])
  }

  first_col = label_cols[1L]
  cluster_assignment[[first_col]] = ifelse(
    cluster_assignment[[first_col]] == "stop", "unknown",
    cluster_assignment[[first_col]]
  )

  cluster_assignment[label_cols] = lapply(
    cluster_assignment[label_cols],
    function(x) ifelse(x == "stop", NA_character_, x)
  )

  cluster_assignment$annotation = apply(
    cluster_assignment[label_cols], 1L,
    function(row) paste(row[!is.na(row)], collapse = "/")
  )

  other_cols = setdiff(names(cluster_assignment), c("cluster_id", "annotation", label_cols))
  cluster_assignment = cluster_assignment[, c("cluster_id", "annotation", label_cols, other_cols)]

  return(cluster_assignment)
}
