# Loading gapminder data and looking at the first couple rows 
gapminder <- read.csv(file="data/gapminder_data.csv")
head(gapminder)

# creating a dataframe from some handwritten data
# c means concatanate, or join into one thing
# so, coat = c("calico", "black", "tabby") concatonates calico, black, and tabby into one object called coat
cats <- data.frame(coat = c("calico", "black", "tabby"),
                   weight = c(2.1, 5.0, 3.2),
                   likes_string = c(1, 0, 1))

# writing some data into a csv file
write.csv(x = cats, file = "data/feline-data.csv", row.names = FALSE)

# Extract the content of the column using $ operator
cats$weight
cats$coat
cats$likes_string

# More operations on columns 
# Add two to every value in a column
cats$weight + 2
# save that as a variable
heavyweight <- cats$weight + 2

# pastes the values of a column into a phrase 
paste("My cat is", cats$coat)
paste("My cat weighs", cats$weight)

# Cannot do a mathematical operation with categorical data and numerical data
#cats$weight + cats$coat

# Command to find out what type of data a column is 
typeof(cats$weight)
typeof(1+1i)
typeof(TRUE)

# The L suffix forces the number to be an integer, since by default R uses float numbers
typeof(1L)


# Another way to view data from a file
file.show(data/felineData.csv)



