#-----------------------------------------------------
#　　繰返し変形による座屈拘束ブレース(BRB)疲労損傷の検証方法
#									CODED by Y.T(2020)
#-----------------------------------------------------
#更新履歴
#2020/05/28		戯れに作る。※塑性ひずみの計算から弾性ひずみは消去しています。
#-----------------------------------------------------
#-----------------------------------------------------
#ModuleのImport
#-----------------------------------------------------
import csv
import os
import sys
#----------------------------------------------------
#Coffin-Manson則の破断繰り返し回数を計算する。
#引数：塑性ひずみ振幅,m2,C2
#戻り値：その振幅の破断繰り返し回数Nf
#----------------------------------------------------
def FuncCoffinMansonNf(dEp,m2,C2):
	Nf = (dEp*100.0 / C2) ** (1.0 / m2)
	return Nf
#----------------------------------------------------
#無次元化軸ひずみ(別名：等価軸ひずみ)時刻歴から極値データを分離する。
#引数：無次元化軸ひずみのリスト(ステップ、歪値)
#戻り値：無次元化軸ひずみの極値リスト(ステップ、歪値)、極値数
#----------------------------------------------------
def FuncKokuchi(E):
	#----------------------------------------------------
	#作業用リストEEにデータの移し替え
	#----------------------------------------------------	
	EE = E
	m = 0
	EEE = []
	dum = 9999.0
	#----------------------------------------------------
	#作業用リストEEの歪が等値な部分はデータを詰める
	#----------------------------------------------------
	for i in EE:
		if i[1] != dum:
			EEE.append(i)
			dum = i[1]
	EE = EEE
	L = len(EE)
	EE[L-1] = EE[-1]
	IP = 0
	Y = [[0.0,0.0]]
	if L < 3:
		print("歪データ数が不足しています。")
		a = input()
		sys.exit()
	#----------------------------------------------------
	#極値の分離
	#----------------------------------------------------		
	for i in range(L-2):
		if EE[i+1][1] > EE[i][1] and EE[i+1][1] > EE[i+2][1]:
			IP = IP + 1
			Y.append(EE[i+1])
		elif EE[i+1][1] < EE[i][1] and EE[i+1][1] < EE[i+2][1]:
			IP = IP + 1
			Y.append(EE[i+1])
	#----------------------------------------------------
	#Return
	#----------------------------------------------------
	return Y,IP
#----------------------------------------------------
#レインフロー法に基づき歪み度分布を計算する。
#引数：軸ひずみ時刻歴の極値リスト,データ数,レインフロー法の範囲数,範囲FR1
#戻り値：頻度数
#----------------------------------------------------
def RainFlow(E,LY,NS,FR):
	IFM = [0.0 for i in range(NS)]
	Y = [line[1] for line in E]
	#----------------------------------------------------
	#データ終端部を極値として処理するための処理
	#----------------------------------------------------
	if Y[-1] != 0.0:
		Y.append(0.0)
		IP = LY + 1
	else:
		IP = LY
	#----------------------------------------------------
	#レインフロー開始
	#----------------------------------------------------
	while IP > 3:
		for i in range(0, IP-3):
			k = 0
			#----増加側----#
			if   Y[i+3] >= Y[i+1] and Y[i+2] >= Y[i]:k=1
			#----減少側----#
			elif Y[i+3] <= Y[i+1] and Y[i+2] <= Y[i]:k=1
			#----レインフロー処理がある場合----#
			if k == 1:
				YY1 = abs(Y[i+3] - Y[i  ])
				YY2 = abs(Y[i+1] - Y[i+2])
				if YY1 >= YY2:
					FM = YY2
					AA = i
					k = 2
			#----変数範囲を決定し間のYを消去して次の極値判定へ----#
			if k == 2:
				for j in range(NS):
					if FM <= FR[j]:
						IFM[j] = IFM[j] + 1.0
						IP = IP - 2
						for l in range(AA,IP-2):
							Y[l] = Y[l+2]
						k=3
						break
				break
		if k <= 1:break
	#----------------------------------------------------
	#小さい順に並べ替える
	#---------------------------------------------------
	Y.sort()
	#---------------------------------------------------
	#両端から振幅を評価する。
	#---------------------------------------------------
	k = IP
	for i in range(int(IP/2)):
		FM = abs(Y[i] - Y[k-i-1])
		for j in range(NS):
			if FM <= FR[j]:
				IFM[j] = IFM[j] + 1.0
				break
	#---------------------------------------------------
	#Return
	#---------------------------------------------------
	return IFM
#----------------------------------------------------
#破断評価のターミナルルーチン
#引数：軸ひずみ時刻歴ファイル名,降伏歪(-),歪増幅係数,as,m2,C,Xs0,損傷度時刻歴の出力フラグ,NS,IS
#戻り値：
#----------------------------------------------------
def FuncFractureBRB(filename,Ey,Amp,As,m2,C2,Xs0,OutFlag,NS,IS):
	#----------------------------------------------------
	#ファイルの有無の確認
	#----------------------------------------------------
	if os.path.exists(filename):
		a = 1
	else:
		print(filename+"が存在しません")
		a = input()
		return
	with open(filename,mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		E = []
		for line in csvline:
			E.append([float(line[0]),float(line[1])*Amp]) #塑性化部を考慮した歪増幅を考慮
	#----------------------------------------------------
	#極値分解
	#----------------------------------------------------
	EE,IP = FuncKokuchi(E)
	#----------------------------------------------------
	#歪み度範囲の作成
	#----------------------------------------------------
	S1 = IS * -0.5			#ある歪度範囲の中央値の計算に使う。
	S2 = (1.0+1.0) * -Ey	#全振幅で見たときは2倍(正+負)の降伏歪を取り除く。
	FR = []
	FR1 = []
	dum = S1 + S2
	Nf = []
	KP = 0.0
	for i in range(NS):
		dum = dum + IS
		FR.append(dum)
		if FR[-1] < 0.0:FR[-1] = 0.0	#レインフロー法で統計した歪度数から弾性分を除いた塑性歪の集計値を出すために使う。
		KP = KP + IS
		FR1.append(round(KP,4))			#レインフロー法で統計する時の歪度範囲は弾性分を含んだ単純な形式
		if FR[i] > 0.0:
			Nf.append(FuncCoffinMansonNf(FR[i],m2,C2*2.0))
		else:
			Nf.append(999999999999.0)
	#----------------------------------------------------
	#フラグの初期化
	#----------------------------------------------------
	YieldFlag = 0
	FractureFlag = 0
	#----------------------------------------------------
	#中間出力ファイルの準備
	#----------------------------------------------------
	if OutFlag == 1:
		file = open('Out_DamageHistory_'+filename,mode='w',newline="")
		out = csv.writer(file)
		out.writerow(['Time(s)','軸ひずみ(無次元)x増幅係数','実効の平均塑性ひずみ(%)','実効の累積塑性ひずみ(%)','疲労破壊が生じる累積塑性ひずみ(%)','判定(1:疲労破壊)','マイナー則で評価した時の損傷度D'])  
	#----------------------------------------------------
	#各ステップの損傷度評価
	#----------------------------------------------------
	CumEp = 0.0	#実効の累積塑性ひずみ
	for i in range(IP):
		if i <= 2: continue
		IFM = RainFlow(EE[0:i+1],i+1,NS,FR1)
		Freq = sum(IFM)
		CumAmp = 0.0	#平均塑性歪を出すための見かけの累積塑性ひずみ
		AveAmp = 0.0	#実効の平均塑性ひずみ
		FratureCumEp = 0.0	#疲労破断判定用の累積塑性ひずみ
		Dminer = 0.0	#マイナー則評価用の損傷度D
		for j in range(NS):CumAmp = CumAmp + FR[j]*IFM[j]	#頻度数に基づく累積歪
		if Freq != 0.0:AveAmp = CumAmp / Freq * 100.0 * 0.5	#平均塑性ひずみ(片振幅)
		if abs(EE[i][1]) >= Ey:YieldFlag = 1				#初期降伏の判定
		if i == 0:
			CumEp = (EE[i][1] - Ey) * 100.0
			if CumEp < 0.0:CumEp = 0.0
		elif YieldFlag == 1:
			CumEp = CumEp + abs(EE[i][1] - EE[i-1][1])*100.0	#実効の累積塑性ひずみ
		FratureCumEp = As / Xs0 + ((1.0-As)*0.25)*((AveAmp**(1.0+m2))/C2)**(-1.0/m2)	#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ
		FratureCumEp = 1.0 / FratureCumEp					#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ
		if FratureCumEp <= CumEp: FractureFlag = 1
		for j in range(NS):
			if FR[j] > 0.0: Dminer = Dminer + IFM[j] / Nf[j]#FuncCoffinMansonNf(FR[j],m2,C2*2.0)
		if OutFlag == 1:
			out.writerow([EE[i][0],EE[i][1],AveAmp,CumEp,FratureCumEp,FractureFlag,Dminer])
	if OutFlag:file.close()
	#----------------------------------------------------
	#return
	#----------------------------------------------------
	return [AveAmp,CumEp,FratureCumEp,FractureFlag,Dminer]

if os.path.exists('InputFileIndex.csv'):
	OutList = []
	with open('InputFileIndex.csv',mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		i = 0
		for i,line in enumerate(csvline):
			if i <= 0:continue
			OutList.append([line[0]]+FuncFractureBRB(line[0],float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),int(line[7]),2000,0.0005))
			print(line[0]+"　評価終了")
	file = open('Out_DamageEvaluationBRB.csv',mode='w',newline="")
	out = csv.writer(file)
	out.writerow(['ファイル名','最終の平均塑性ひずみ片振幅(%)','最終の累積塑性ひずみ(%)','疲労破壊が生じる累積塑性ひずみ(%)','判定(1:疲労破壊)','マイナー則で評価した時の損傷度D'])
	for line in OutList:
		out.writerow(line)
	file.close()
else:
	print('InputFileIndex.csv'+"が存在しません")
	a=input()