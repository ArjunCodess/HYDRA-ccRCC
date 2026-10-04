# Verify the checkpoint writer preserves the exact stage-22 null mechanism.
get_assignment <- function(expressions, variable) {
  found <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    is.symbol(x[[2]]) && identical(x[[2]], as.name(variable)), as.list(expressions))
  stopifnot(length(found) == 1L)
  found[[1]]
}
primary <- parse("analysis/22_nested_cv.R", keep.source = FALSE)
checkpointed <- parse("analysis/50_clinical_null_refits.R", keep.source = FALSE)
for (name in c("baseline", "null_lp", "base_hazard", "censor_times")) {
  stopifnot(identical(get_assignment(primary, name), get_assignment(checkpointed, name)))
}
find_loop <- function(x) {
  if (is.call(x) && identical(x[[1]], as.name("for")) && identical(x[[2]], as.name("sim"))) return(list(x))
  if (!is.recursive(x)) return(list())
  unlist(lapply(as.list(x), find_loop), recursive = FALSE)
}
original_loop <- find_loop(primary)[[1]]
new_loop <- find_loop(checkpointed)[[1]]
for (name in c("simulated_event_time", "simulated_censor", "sim_time", "sim_event", "fold", "f", "selected", "sim_pred")) {
  stopifnot(identical(get_assignment(as.list(original_loop[[4]])[-1], name),
                      get_assignment(as.list(new_loop[[4]])[-1], name)))
}
# Restoring the captured RNG state must reproduce later draws, rather than
# restarting the seed and generating duplicate simulated patients/outcomes.
set.seed(20260749)
expected <- replicate(20, c(runif(10), sample(1:12, 10, replace = TRUE)))
set.seed(20260749)
first <- replicate(7, c(runif(10), sample(1:12, 10, replace = TRUE)))
saved_rng <- .Random.seed
runif(17)
assign(".Random.seed", saved_rng, envir = .GlobalEnv)
remaining <- replicate(13, c(runif(10), sample(1:12, 10, replace = TRUE)))
stopifnot(identical(cbind(first, remaining), expected))
message("PASS: original clinical-null DGP and selector calls unchanged; split/resumed RNG stream matches uninterrupted draws.")
