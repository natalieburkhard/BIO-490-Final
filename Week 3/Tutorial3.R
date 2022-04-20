# Loading gapminder data and looking at the first couple rows 
gapminder <- read.csv(file="data/gapminder_data.csv")
head(gapminder)

# loading the package ggplot2 for making graphs
library("ggplot2")
# plotting a scatter plot. arguments of the function ggplot are (data, mapping) data tells it what data to use and mapping tells it what it will look like and what variable goes on what axis
# think of the geom_point() as putting the datapoints on the graph
ggplot(data = gapminder, mapping = aes(x = gdpPercap, y = lifeExp)) +
  geom_point()

#plotting life expectancy vs time
ggplot(data = gapminder, mapping = aes(x = year, y = lifeExp)) + geom_point()

# changing the colors of the data points by continent
ggplot(data = gapminder, mapping = aes(x = year, y = lifeExp, color = continent)) + geom_point()

# plotting a line plot
# using geom_line() instead of geo_point to make the line plot
ggplot(data = gapminder, mapping = aes(x=year, y=lifeExp, color=continent)) +
  geom_line()

# plotting each country seperately (adding the argument by)
ggplot(data = gapminder, mapping = aes(x=year, y=lifeExp, by=country, color=continent)) +
  geom_line()

# adding datapoints to our line graph (we're adding another layer here)
ggplot(data = gapminder, mapping = aes(x=year, y=lifeExp, by=country, color=continent)) +
  geom_line() + geom_point()

ggplot(data = gapminder, mapping = aes(x=year, y=lifeExp, by=country)) +
  geom_line(mapping = aes(color=continent)) + geom_point()

# switch the order of the point and line layers 
ggplot(data = gapminder, mapping = aes(x=year, y=lifeExp, by=country)) +
   geom_point() + geom_line(mapping = aes(color=continent))

# transforming to a logarithmic scale
ggplot(data = gapminder, mapping = aes(x = gdpPercap, y = lifeExp)) +
  geom_point(alpha = 0.5) + scale_x_log10()

# fitting a line to the plot 
ggplot(data = gapminder, mapping = aes(x = gdpPercap, y = lifeExp)) +
  geom_point(alpha = 0.5) + scale_x_log10() + geom_smooth(method="lm")

# changing the size of the line
ggplot(data = gapminder, mapping = aes(x = gdpPercap, y = lifeExp)) +
  geom_point(alpha = 0.5) + scale_x_log10() + geom_smooth(method="lm", size=1.5)

# changing it so the points are colored by continent and change the shape
ggplot(data = gapminder, mapping = aes(x = gdpPercap, y = lifeExp, color = continent)) +
  geom_point(size=3, shape=17) + scale_x_log10() +
  geom_smooth(method="lm", size=1.5)

# Multipanel figures
#americas <- gapminder[gapminder$continent == "Americas",]
#ggplot(data = americas, mapping = aes(x = year, y = lifeExp)) +
#  geom_line() +
#  facet_wrap( ~ country) +
#  theme(axis.text.x = element_text(angle = 45))

