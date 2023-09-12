test_that("Correct number of columns in output", {
  expect_equal(dim(AssignmentAgreement(SDT))[2], 18)
})