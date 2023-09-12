test_that("error if rarefy depth too deep", {
  expect_error(GetMultiRarefiedRichnessEstimates(ps=ps, RarefyDepth=50000, NumReplicates=2),)
})

test_that("correct dimensions of output", {
  expect_equal(dim(GetMultiRarefiedRichnessEstimates(ps=ps, RarefyDepth=100, NumReplicates=3)), c(6,7))
})

