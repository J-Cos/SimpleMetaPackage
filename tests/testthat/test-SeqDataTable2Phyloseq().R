SDTfile <- system.file("extdata", "16s_multirun_test_SeqDataTable.RDS", package="SimpleMetaPackage")
deprecatedSDTfile <- system.file("extdata", "COI_test_SeqDataTable.RDS", package="SimpleMetaPackage")


test_that("produces phyloseq", {
  expect_equal(class(SeqDataTable2Phyloseq(SDTfile))[1], "phyloseq")
})

test_that("all clustering options work", {
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "ESV")), 110)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedESV")), 101)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "OTU")), 97)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU")), 47)
    expect_error(SeqDataTable2Phyloseq(SDTfile, clustering = "badOption"))
})

test_that("all clustering options work", {
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "ESV")), 110)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedESV")), 101)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "OTU")), 97)
    expect_equal(ntaxa(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU")), 47)
})

test_that("different assignment proportions gets different results", {
    expect_false( isTRUE( all.equal( tax_table(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment="RepresentativeSequence")), tax_table(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment=0.1)) )))
    expect_false( isTRUE( all.equal( tax_table(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment=0.8)), tax_table(SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment=0.1)) )))
})

test_that("assignment proportions >1 fail", {
    expect_error( SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment=1000) )
    expect_error( SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment=-0.1) )
    expect_error( SeqDataTable2Phyloseq(SDTfile, clustering = "curatedOTU", ClusterAssignment="other") )
})


test_that("standard Fastq naming argument works", {
    expect_identical( sort(sample_names(SeqDataTable2Phyloseq(SDTfile, StandardFastqNaming=FALSE))),  c("SMALLDuplicate_RVS_Low_Illi_123_10_L001", "SMALL_CS_Low_Illi_127_42_L001", "SMALL_RVS_Low_Illi_123_10_L001"))
    expect_identical( sort(sample_names(SeqDataTable2Phyloseq(SDTfile, StandardFastqNaming=TRUE))),  c("SMALLDuplicate_RVS_Low_Illi_123", "SMALL_CS_Low_Illi_127", "SMALL_RVS_Low_Illi_123"))
})

test_that("expect warning on deprecated SDTs (before multi run management implemented) if non standard fastq naming chosen", {
    expect_warning( SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=FALSE))
})

test_that("expect error on deprecated SDTs (before multi run management implemented) if standard fastq naming chosen ", {
    expect_error( SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=TRUE))
})

test_that("BLAST options work", {
    expect_false( isTRUE( all.equal( SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=FALSE, assignment="BLAST"), SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=FALSE, assignment="Idtaxa") )))
    expect_false( isTRUE( all.equal(  SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=FALSE, assignment="BLAST", BLASTThreshold=80),  SeqDataTable2Phyloseq(deprecatedSDTfile, StandardFastqNaming=FALSE, assignment="BLAST", BLASTThreshold=1) )))
})

test_that("no assignment works", {
  expect_equal(class(SeqDataTable2Phyloseq(SDTfile, assignment="None"))[1], "phyloseq")
  expect_equal(class(SeqDataTable2Phyloseq(deprecatedSDTfile, assignment="None", StandardFastqNaming=FALSE))[1], "phyloseq")
})

test_that("invalid assignment throws error", {
  expect_error(SeqDataTable2Phyloseq(SDTfile, assignment="badOption"))
})