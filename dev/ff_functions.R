get_tc_vars <- function() {
  return(c('c_sex_female','c_lbw','f_age','f_immigrant','m_age','m_immigrant','m_living_parents_15','m_cognitive'))
}

get_tv_vars <- function() {
  return(c('m_f_in_jail','m_kids','m_evicted','m_anxiety','m_depression','m_health_excellent','m_health_very_good','m_health_good',
           'm_health_fair','m_health_poor','m_job_type_full_time','m_job_type_part_time','m_job_type_unemployed',
           'm_relation_married','m_relation_cohab','m_relation_contact','m_relation_no_contact','m_wodtke'))
}

clean_cov_names <- function() {
  
  ## Variable
  tc_vars <- get_tc_vars()
  tv_vars <- get_tv_vars()
  dvs <- c('c_ppvt')
  
  ## Clean names
  tc_var_names <- c('Percent female','Low birth weight','Father age','Father immigrant status','Mother age','Mother immigrant status',
                    'Mother living with parents at age 15','Mother cognitive score')
  tv_var_names <- c('Father incarcerated','Total children','Mother evicted','Mother anxiety','Mother depression','Excellent','Very good',
                    'Good','Fair','Poor','Full time','Part time','Unemployed','Married','Cohabitating','Contact','No contact',
                    'Wodtke index')
  dv_names <- c('Child PPVT score')
  
  dt <- data.table(variable = c(tc_vars, tv_vars, dvs, 'censor'),
                   name = c(tc_var_names, tv_var_names, dv_names, 'Censoring'))
  dt[, cov_sort := 1:.N]
  return(dt)
  
}