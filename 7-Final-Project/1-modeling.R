options(scipen = 999)
glm_mod_1 <- glm(cbind(tri_np, round(tri_nneg)) ~
                   alt + bio12,
                 data=STH_kenya_df, family=binomial())
glm_mod_1 %>% summary()


glm_mod_2 <- glm(cbind(tri_np, round(tri_nneg)) ~
                   alt + bio12 + lakes,
                 data=STH_kenya_df, family=binomial())
glm_mod_2 %>% summary()