#run_ard_compare_grouped_ieugwasr
system("git fetch origin")
system("git pull origin main")

setwd('G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\ardmr')

devtools::load_all()
rm(list = ls())

#set token
#Sys.setenv(IEU_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJoLnJvYmVydC5jYW1wYmVsbEBnbWFpbC5jb20iLCJpYXQiOjE3NTg1ODE1NjksImV4cCI6MTc1OTc5MTE2OX0.DoSCsFOVcrZZtQp9uOrURmnXzjX8Zm6iVwVMAwuzvW5rVgKkepKgKyasJDHa4eDi7IgEJqH0thBzI30IOW4__mj59RkXyBAFXfoQ-Gfhnjvln8A54gTlhc8671MIYmvOba8491nIHYgWe1PgWYv9JXJ_9X35YQ2GNtPa9i49DK-t3yyhwBJ3UrIPTiOy8zYx_gzJOB9pcTb7hs6AMC9UotcqzHvJVlRpED3ADYOq2UQE3HsbzR3hbGo3kzfLLPXashtjbni_sHZZW1qZsyKNNN3qbQgnbV3-ju5_R9ElYg5UNaq3AbNlcoBNsNxkbv5ojPjVaCdNkBVANV7mNqbwbA")
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJoLnJvYmVydC5jYW1wYmVsbEBnbWFpbC5jb20iLCJpYXQiOjE3NTg1ODE1NjksImV4cCI6MTc1OTc5MTE2OX0.DoSCsFOVcrZZtQp9uOrURmnXzjX8Zm6iVwVMAwuzvW5rVgKkepKgKyasJDHa4eDi7IgEJqH0thBzI30IOW4__mj59RkXyBAFXfoQ-Gfhnjvln8A54gTlhc8671MIYmvOba8491nIHYgWe1PgWYv9JXJ_9X35YQ2GNtPa9i49DK-t3yyhwBJ3UrIPTiOy8zYx_gzJOB9pcTb7hs6AMC9UotcqzHvJVlRpED3ADYOq2UQE3HsbzR3hbGo3kzfLLPXashtjbni_sHZZW1qZsyKNNN3qbQgnbV3-ju5_R9ElYg5UNaq3AbNlcoBNsNxkbv5ojPjVaCdNkBVANV7mNqbwbA")

#set cache dir
Sys.setenv(ARDMR_CACHE_DIR = 'C:\\Users\\Robert\\Downloads\\ardmr')

#set csv path
csv_path = 'G:\\My Drive\\Documents\\0Oxford_main\\ARD paper\\3_ARD_MR\\exposures to run on\\ieugwasr_exposures_23.9.25_ancestry_specific_only.csv'

#call the ard_compare_grouped_ieugwasr
result <- run_ieugwasr_ard_compare(
  csv_path = csv_path,
  cache_dir = ardmr_cache_dir(),   # or your custom path
  prompt_for_units = TRUE,          # will ask you for any missing units
  sensitivity_pass_min = 5,
  force_refresh = TRUE
  )
