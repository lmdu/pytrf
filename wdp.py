seq = "AAGCCGAGAAGGTAGATAG"
motif = "AAG"
dp = []

n = len(seq)
m = len(motif)

for i in range(n+1):
	dp.append([])
	for j in range(m+1):
		if i == 0:
			dp[i].append(j)
		elif j == 0:
			dp[i].append(i)
		else:
			dp[i].append(0)

for i in range(1, n+1):
	#first pass
	cost = 0 if motif[0] == seq[i-1] else 1
	dp[i][1] = min(dp[i-1][0]+cost,dp[i-1][m]+cost, dp[i-1][1]+1)
	for j in range(2, m+1):
		cost = 0 if motif[j-1] == seq[i-1] else 1
		dp[i][j] = min(dp[i-1][j-1]+cost, dp[i][j-1]+1, dp[i-1][j]+1)

	#second pass
	dp[i][1] = min(dp[i][1], dp[i][m]+1)
	for j in range(2, m):
		dp[i][j] = min(dp[i][j], dp[i][j-1]+1)

for row in dp:
	print(*row, sep='\t')

