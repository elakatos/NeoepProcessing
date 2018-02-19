setwd('~/Software/NeoepProcessing/')

test_that('avinput files are read and named correctly',
          {
          exampleFile <- 'test/example1.avinput'
          avinput <- readAvinput(exampleFile)
          expect_equal('Region2', names(avinput)[20])
            })
