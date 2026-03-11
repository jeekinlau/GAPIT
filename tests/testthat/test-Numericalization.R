test_that("GAPIT.Numericalization.Batch matches serial reference across scenarios", {
  bit2_x <- rbind(
    c("AA", "AT", "TT", "NN", "AA", "AT"),
    c("CC", "CT", "TT", "NN", "CT", "CC"),
    c("GG", "GG", "GG", "NN", "GG", "NN"),
    c("AT", "TA", "AA", "TT", "NN", "AA")
  )

  bit1_x <- rbind(
    c("A", "R", "T", "N", "A", "R"),
    c("C", "Y", "T", "N", "Y", "C"),
    c("G", "G", "G", "N", "G", "N"),
    c("R", "A", "T", "C", "N", "A")
  )

  scenarios <- expand.grid(
    bit = c(1, 2),
    effect = c("Add", "Hybrid", "Dominant", "Recessive"),
    impute = c("Middle", "Minor", "Major"),
    Major.allele.zero = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    x <- if (s$bit == 1) bit1_x else bit2_x

    expected <- apply(x, 1, function(one) {
      GAPIT.Numericalization(
        x = one,
        bit = s$bit,
        effect = s$effect,
        impute = s$impute,
        Major.allele.zero = s$Major.allele.zero
      )
    })

    actual_serial <- GAPIT.Numericalization.Batch(
      X = x,
      bit = s$bit,
      effect = s$effect,
      impute = s$impute,
      Major.allele.zero = s$Major.allele.zero,
      ncpus = 1
    )

    actual_parallel <- GAPIT.Numericalization.Batch(
      X = x,
      bit = s$bit,
      effect = s$effect,
      impute = s$impute,
      Major.allele.zero = s$Major.allele.zero,
      ncpus = 2
    )

    expect_equal(actual_serial, expected)
    expect_equal(actual_parallel, expected)
  }
})


test_that("GAPIT.HapMap serial and parallel numericalization are identical", {
  g <- as.data.frame(
    rbind(
      c("rs#", "alleles", "chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panelLSID", "QCcode", "taxa1", "taxa2", "taxa3"),
      c("snp1", "A/T", "1", "100", "+", "NA", "NA", "NA", "NA", "NA", "NA", "AA", "AT", "TT"),
      c("snp2", "C/T", "1", "200", "+", "NA", "NA", "NA", "NA", "NA", "NA", "CC", "CT", "NN"),
      c("snp3", "G/T", "2", "300", "+", "NA", "NA", "NA", "NA", "NA", "NA", "GG", "GG", "GT")
    ),
    stringsAsFactors = FALSE
  )

  hm_serial <- GAPIT.HapMap(g, SNP.effect = "Add", SNP.impute = "Middle", heading = TRUE, ncpus = 1)
  hm_parallel <- GAPIT.HapMap(g, SNP.effect = "Add", SNP.impute = "Middle", heading = TRUE, ncpus = 2)

  expect_equal(hm_parallel$GD, hm_serial$GD)
  expect_equal(hm_parallel$GT, hm_serial$GT)
  expect_equal(hm_parallel$GI, hm_serial$GI)
})


test_that("GAPIT.HapMap indicator mode remains unchanged", {
  g <- as.data.frame(
    rbind(
      c("snp1", "A/T", "1", "100", "+", "NA", "NA", "NA", "NA", "NA", "NA", "AA", "AT", "TT"),
      c("snp2", "C/T", "1", "200", "+", "NA", "NA", "NA", "NA", "NA", "NA", "CC", "CT", "NN")
    ),
    stringsAsFactors = FALSE
  )

  hm_indicator <- GAPIT.HapMap(g, heading = FALSE, Create.indicator = TRUE, ncpus = 2)
  expect_equal(hm_indicator$GD, t(g[, -(1:11)]))
})
