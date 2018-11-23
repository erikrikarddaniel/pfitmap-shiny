app <- ShinyDriver$new("../")
app$snapshotInit("mytest")

# Input 'mainmatrix_rows_current' was set, but doesn't have an input binding.
# Input 'mainmatrix_rows_all' was set, but doesn't have an input binding.
app$snapshot()
