require(compareGroups)

train = readRDS('~/git/IBDpred/01_Preprocessing/data/train__IBD.rds')
test = readRDS('~/git/IBDpred/01_Preprocessing/data/test__IBD.rds')


require(phyloseq)
get_sample(train, 'IBD')

train.df = data.frame(
  cohort = 'train',
  Age = get_variable(train, 'AGE_YEARS'),
  Sex = get_variable(train, 'SEX'),
  IBD = ifelse(get_variable(train, 'IBD') == 'YES', 'IBD', 'control')
)

test.df = data.frame(
  cohort = 'test',
  Age = get_variable(test, 'AGE_YEARS'),
  Sex = get_variable(test, 'SEX'),
  IBD = ifelse(get_variable(test, 'IBD') == 'YES', 'IBD', 'control')
)

data = rbind.data.frame(train.df, test.df)

require(compareGroups)
cg = compareGroups(cohort ~., data = data)
res = createTable(cg)

export2latex(res)
require(xtable)
xtable(res)

