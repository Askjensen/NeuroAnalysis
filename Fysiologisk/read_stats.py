import pstats
p = pstats.Stats('profile.result')
p.strip_dirs().sort_stats(0).print_stats()