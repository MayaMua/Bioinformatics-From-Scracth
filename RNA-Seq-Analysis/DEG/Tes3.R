install.packages("gapminder")

library(gapminder)

countries <- gapminder
str(countries)
install.packages("dplyr")
install.packages("ggplot2")

library(dplyr)
library(ggplot2)
countries %>% 
  filter(year == 1952)
countries %>% 
  filter(country == "Pakistan")
countries %>%
  arrange(lifeExp)
countries %>% 
  arrange(desc(lifeExp))
countries %>%
  mutate(pop = pop / 1000000)

countries %>% 
  filter(country == "Pakistan") %>%
  arrange(desc(lifeExp)) %>%
  mutate(pop = pop / 1000000)

ggplot(countries, aes(x=gdpPercap, y=lifeExp)) + geom_point()

