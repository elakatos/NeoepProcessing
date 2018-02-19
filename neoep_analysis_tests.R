setwd('~/Software/NeoepProcessing/')

test_that('avinput files are read and named correctly',
          {
          exampleFile <- 'test/example1.avinput'
          avinput <- readAvinput(exampleFile)
          expect_equal('Region2', names(avinput)[20])
            })

test_that('region names are generated correctly',
          {
            expect_equal(c('Region1', 'Region2'), getRegionNames(2))
            expect_equal(c('Normal','Region1','Region2','Region3','Region4','Region5'), getRegionNames(5, T))
          })

test_that('vaf computation is correct',
          {
            data <- data.frame(c("10:0", "10:1", "10:2", "20:2", "20:10"), stringsAsFactors = F)
            expect_equal(c(0, 0.1, 0.2, 0.1, 0.5), computeVaf(data, 1))
          })
