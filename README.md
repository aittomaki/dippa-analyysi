dippa-analyysi
==============
Repo dippatyön analyysejä varten. Tänne tulee kaikki koodit yms.


Repon sisältö
-------------

- README.md -- tämä helppitiedosto.
- dippa.and -- pääanalyysiscripti
- dippa-analyysi.Rproj -- Rstudio-projekti datalla ja tuloksilla leikkimiseen
- preprocess.R -- DEPRECATED! PREPROSESSOINTI ANDURILILLA. R-scriptin datan esikäsittelyyn.
- regression_analysis_multivariate.R -- R-scripti, jossa Bayesilainen regressioanalyysi
  (prot ~ mrna + mirna), tässä kaikki miRNAt kerralla mukana.
- regression_analysis_univariate.R -- R-scripti, jossa Bayesilainen regressioanalyysi
  (prot ~ mrna + mirna), tässä yksi miRNA kerrallaan mallissa.
- scale_regression_data.R -- DEPRECATED! NORMALISOINTI ANDURILILLAR-scripti datan 
  normalisointiin regressiota varten. Plottaa myös jakaumat muuttujista.
- shrinkage_prior.stan -- Bayes-mallitiedosto Stanilla simulointia varten, tässä
  mallissa HS-priori miRNAlle.
- simple_priors.stan -- Bayes-mallitiedosto Stanilla simulointia varten. Tässä
  tiedostossa yksinkertaisin mahdollinen malli (erityisesti priorit)-
- testing_data.R -- R-scipti, jossa lähinnä datan kanssa leikkimistä ja uuden
  (norjalaisten paperissa julkaistun) ja vanhan (norjalaisten mulle lähettämän)
  datan vertaamista
