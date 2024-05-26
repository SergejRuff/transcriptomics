# Generate the matrix
mat <- matrix(c(0, 2, 0, 1, 2,
                1, 3, 0, 1, 0,
                1, 1, 1, 3, 1,
                2, 1, 1, 1, 0,
                0, 0, 2, 2, 0), nrow = 5, byrow = TRUE)

mat <- as(mat, "sparseMatrix")

# Print the matrix
print(mat)
# Define a scale factor
scale.factor <- 1



# Compute column pointers manually
col_pointers <- c(1, cumsum(Matrix::colSums(mat)))

# Compute differences between consecutive pointers
pointer_diff <- diff(col_pointers)

# Repeat the column sums according to pointer_diff
repeated_col_sums <- rep.int(Matrix::colSums(mat), pointer_diff)

# Perform normalization
norm.data <- mat
norm.data <- norm.data / repeated_col_sums * scale.factor

# Display the normalized matrix
print("Normalized Matrix:")
print(norm.data)

## selbes Ergebnis wie: https://satijalab.org/seurat/reference/relativecounts. Also richtig