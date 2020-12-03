# Ascaris lumbricoides (roundworm)
STH_kenya %>% ggplot(aes(x=asc_prev)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666")

STH_kenya %>% ggplot(aes(x=asc_prev_logit)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666")

# Hookworm
STH_kenya %>% ggplot(aes(x=hk_prev)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666") 

STH_kenya %>% ggplot(aes(x=hk_prev_logit)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666") 

# roundworm
pal = colorNumeric("Oranges", STH_kenya$asc_prev_logit) 

STH_kenya %>%
  leaflet() %>% addProviderTiles("CartoDB.DarkMatter") %>%
  addCircleMarkers(~long, ~lat, fillOpacity=1,
                   fillColor= ~pal(asc_prev_logit),
                   radius=3,
                   weight=0.1,
                   stroke=TRUE) %>%
  addLegend(pal = pal, values = ~asc_prev_logit, title = "Roundworm Prevalence (%)")

# hookworm
pal = colorNumeric("Oranges", STH_kenya$hk_prev_logit) 

STH_kenya %>%
  leaflet() %>% addProviderTiles("CartoDB.DarkMatter") %>%
  addCircleMarkers(~long, ~lat, fillOpacity=1,
                   fillColor= ~pal(hk_prev_logit),
                   radius=3,
                   weight=0.1,
                   stroke=TRUE) %>%
  addLegend(pal = pal, values = ~asc_prev_logit, title = "Hookworm Prevalence (%)")
