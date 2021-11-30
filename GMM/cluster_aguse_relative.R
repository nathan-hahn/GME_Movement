df <- output

table(df$ag.use, df$tactic.season)

prop.table(table(df$ag.use, df$tactic.season), 1)
