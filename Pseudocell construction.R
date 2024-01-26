AvgN <- function(data, cell.info, N=10, seed=1024, field="typeCourse") {
  cell.types <- table(cell.info[[field]])
  ## drop the cell states less than 20 cells
  cell.types <- cell.types[cell.types >= N]
  cell.types <- sort(cell.types, decreasing = TRUE)
  new.data <- lapply(names(cell.types), function(cc) {
    ## get cells in same cell types and tissue (cell state)
    new.cell.info <- subset(cell.info, get(field) == cc)
    ## set the seed to make the shuffle function reproducible
    set.seed(seed)
    new.cell.info <- new.cell.info[shuffle(rownames(new.cell.info)), ]
    ## generate PoolID for pooling
    new.cell.info$PoolID <- floor(0:(nrow(new.cell.info)-1) / N)
    ## pooling cells with the same PoolID
    pool.ids <- table(new.cell.info$PoolID)
    pool.ids <- pool.ids[pool.ids >= N]
    pool.ids <- sort(names(pool.ids))
    ## average the expression profiles
    avgN.data <- lapply(pool.ids, function(x) {
      select.cells <- rownames(subset(new.cell.info, PoolID == x))
      new.data <- data[, select.cells]
      return(rowMeans(new.data))
    })
    avgN.data <- do.call(cbind, avgN.data)
    colnames(avgN.data) <- paste(cc, pool.ids, sep = ".")
    return(avgN.data)
  })
  ## merge the expression profiles for each cell state
  new.data <- do.call(cbind, new.data)
  return(new.data)
}