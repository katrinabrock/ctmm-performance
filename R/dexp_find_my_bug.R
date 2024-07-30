
dexp2_old <- function(x,Exp=exp(-x)) { ifelse(Exp<0.7071068,1-Exp^2,2*Exp*sinh(x)) }
# 1 - exp(-x)^1
dexp1_old <- function(x,Exp=exp(-x)) { ifelse(Exp<0.5,1-Exp,2*sqrt(Exp)*sinh(x/2)) }

dexp2_new <- function(x,Exp=exp(-x)) { if(Exp<0.7071068) return(1-Exp^2) else return(2*Exp*sinh(x))}
# 1 - exp(-x)^1
dexp1_new <- function(x,Exp=exp(-x)) { if(Exp<0.1) return(1-Exp) else return(2*sqrt(Exp)*sinh(x/2))}




test_cases <- list(
                   


structure(c(2, 0.495140915152137, 0.609485015571034), names = c("", "", "")),
structure(c(2, 0.495021136306656, 0.609558023354844), names = c("", "", "")),
structure(c(2, 0.453124999910078, 0.63563867388321), names = c("", "", "")),
c(2, position = 0.00648069023167319, position = 0.993540264150544),
structure(c(2, 0.495077071427664, 0.609523928606602), names = c("", "", "")),
c(2, position = 0.00694529209767019, position = 0.993078770703536),
structure(c(2, 0.033565082407907, 0.966991975017244), names = c("", "", "")),
structure(c(2, 0.52025048183328, 0.59437165002231), names = c("", "", "")),
structure(c(2, 0.520254629542668, 0.59436918474655), names = c("", "", "")),
structure(c(2, 0.528581891302551, 0.589440267636884), names = c("", "", "")),
c(2, position = 0.00717504989917924, position = 0.992850629318089),
c(2, position = 0.00729166649779725, position = 0.992734853205602),
c(2, position = 0.00694359667667722, position = 0.993080454391559),
structure(c(2, 0.511863425871563, 0.599377641734874), names = c("", "", "")),
structure(c(2, 0.0335688755955864, 0.966988307042155), names = c("", "", "")),
c(2, position = 0.0067129628076242, position = 0.993309518793025),
c(2, position = 0.00717680183420529, position = 0.992848889909819),
structure(c(2, 0.486685934625792, 0.614660046071054), names = c("", "", "")),
c(2, position = 0.000462962952542656, position = 0.999537144198269),
structure(c(2, 0.511863425842257, 0.599377641752439), names = c("", "", "")),
structure(c(2, 0.469907407340205, 0.625060141584379), names = c("", "", "")),
structure(c(2, 0.50341132505563, 0.60446511160724), names = c("", "", "")),
structure(c(2, 0.478302424429899, 0.61983471547347), names = c("", "", "")),
c(2, position = 0.00648148133115076, position = 0.993539478161671),
c(2, position = 0.0138871933533544, position = 0.986208788894545),
c(2, position = 0.00671378236113508, position = 0.993308704723055),
c(2, position = 0.00659802749286753, position = 0.993423691696342),
c(2, position = 0.00648148133155498, position = 0.993539478161269),
structure(c(2, 0.495084965633853, 0.609519116918025), names = c("", "", "")),
structure(c(2, 0.461516203639795, 0.63032721619989), names = c("", "", "")),
structure(c(2, 0.469911153827922, 0.625057799808623), names = c("", "", "")),
c(2, position = 0.000347222214103827, position = 0.999652838060553),
structure(c(2, 0.0335645472074875, 0.966992492551893), names = c("", "", "")),
structure(c(2, 0.0335648148309475, 0.966992233762051), names = c("", "", "")),
c(2, position = 0.0138905841953404, position = 0.986205444822046),
c(2, position = 0.0065964168429242, position = 0.993425291755443),
structure(c(2, 0.0419560185020519, 0.95891195402919), names = c("", "", "")),
structure(c(2, 0.520254629571973, 0.594369184729132), names = c("", "", "")),
structure(c(2, 0.486689814741026, 0.614657661123873), names = c("", "", "")),
c(2, position = 0.00705932328792823, position = 0.992965535205604),
structure(c(2, 0.48668981471172, 0.614657661141886), names = c("", "", "")),
c(2, position = 0.00682953722940264, position = 0.9931937310593),
c(2, position = 0.00717592576016686, position = 0.992849759719337),
structure(c(2, 0.528645833243078, 0.589402578887305), names = c("", "", "")),
structure(c(2, 0.48674869625324, 0.61462147021679), names = c("", "", "")),
c(2, position = 0.0138905797089655, position = 0.986205449246533),
c(2, velocity = 1.99998225320918, velocity = 0.135337685024887),
c(1, 1.00693641645907, position = 0.133470762459955)
)

for(test_case in test_cases){
    if (
        list(
            dexp1_old,
            dexp2_old
        )[[test_case[1]]](test_case[2], test_case[3]) != list(
            dexp1_new,
            dexp2_new
        )[[test_case[1]]](test_case[2], test_case[3])
        ) dput(test_case)
}

corner_case <- c(1, 1.00693641645907, position = 0.133470762459955)
dexp1_old(corner_case[2], corner_case[3])
dexp1_new(corner_case[2], corner_case[3])
