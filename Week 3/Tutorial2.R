
# creating a simple numeric vector:
x <- c(5.4, 6.2, 7.1, 4.8, 7.5)
names(x) <- c('a', 'b', 'c', 'd', 'e')
x

# extract elements of a vector by giving their corresponding index, starting from one
# the square brackets operator is a function. For vectors (and matrices), it means “get me the nth element”.
x[1]

# retrieve multiple elemnts at the same time
x[c(1, 3)]

# retrieving a slice of the vector
# the : operator creates a sequence of numbers from the left element to the right.
x[1:4]
1:4
c(1, 2, 3, 4)

# retireve the same element multiple times
x[c(1,1,3)]

# can't ask for an index outside the length of the vector. It gives us a vector of length one containing an NA, whose name is also NA.
x[6]

# trying to retrieve the 0th element returns an empty vector
x[0]
# Don't get tripped up! in Python and C the first element of a vector has an index of 0. In R, the first element is 1.

# If we use a negative number as the index of a vector, R will return every element except for the one specified:
x[-2]
x[c(-1, -5)]

# : is really a function. It takes its first argument as -1, and its second as 3, so generates the sequence of numbers: c(-1, 0, 1, 2, 3)
# so you can't skip slices of a vector
# you have to wrap that function call in brackets, so that the - operator applies to the result
x[-(1:3)]
