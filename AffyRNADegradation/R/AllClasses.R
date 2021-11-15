setClass("AffyDegradationBatch", representation(location.type = "character",
                                                afbatch       = "AffyBatch",
                                                stats         = "data.frame",
                                                means.pm      = "matrix",
                                                means.mm      = "matrix"))
