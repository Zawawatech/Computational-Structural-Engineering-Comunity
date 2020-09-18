#-----------------------------------------------------
#　　繰返し変形による座屈拘束ブレース(BRB)疲労損傷の検証方法
#									CODED by Y.T(2020)
#-----------------------------------------------------
#更新履歴
#2020/05/28		戯れに作る。※塑性ひずみの計算から弾性ひずみは消去しています。
#-----------------------------------------------------
#Copyright (c) <2020>, <Yuki TERAZAWA>
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met: 
#
#1. Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer. 
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution. 
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#The views and conclusions contained in the software and documentation are those
#of the authors and should not be interpreted as representing official policies, 
#either expressed or implied, of the FreeBSD Project.
#-----------------------------------------------------
#ModuleのImport
#-----------------------------------------------------
import csv
import os
import sys
#----------------------------------------------------
#Coffin-Manson則の破断繰り返し回数を計算する。
#引数：塑性ひずみ振幅(無次元),m2,C2
#戻り値：その振幅の破断繰り返し回数Nf
#----------------------------------------------------
def FuncCoffinMansonNf(dEp,m2,C2):
	Nf = (dEp*100.0 / C2) ** (1.0 / m2)
	return Nf
#----------------------------------------------------
#時刻歴から極値データを分離する。
#引数：入力データままのリスト(ステップ、歪値)
#戻り値：極値リスト(ステップ、歪値)、極値数、経験最大歪のリスト
#----------------------------------------------------
def FuncKokuchi(E):
	#----------------------------------------------------
	#作業用リストEEにデータの移し替え
	#----------------------------------------------------	
	EE = E
	m = 0
	EEE = []
	dum = 9999.0
	AS = []
	Emax = 0.0
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
	AS = [0.0]
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
			if abs(EE[i+1][1]) > Emax:Emax = abs(EE[i+1][1])
			AS.append(Emax*100.0)
			if i == L-3:
				IP = IP + 1
				Y.append(EE[i+2])
				if abs(EE[i+2][1]) > Emax:Emax = abs(EE[i+2][1])
				AS.append(Emax*100.0)
		elif EE[i+1][1] < EE[i][1] and EE[i+1][1] < EE[i+2][1]:
			IP = IP + 1
			Y.append(EE[i+1])
			if abs(EE[i+1][1]) > Emax:Emax = abs(EE[i+1][1])
			AS.append(Emax*100.0)
			if i == L-3:
				IP = IP + 1
				Y.append(EE[i+2])
				if abs(EE[i+2][1]) > Emax:Emax = abs(EE[i+2][1])
				AS.append(Emax*100.0)
	#----------------------------------------------------
	#Return
	#----------------------------------------------------
	return Y,IP,AS
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
	#---------------------------------------------------
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
#引数：(ファイル名,増幅フラグ,降伏歪(-),塑性下部長さ比,塑性化部断面積比,m2,C2,OutFlag)
#戻り値：評価結果
#----------------------------------------------------
def FuncFractureBRB(filename,AmpFlag,Ey,LPL0,APAE,m2,C2,OutFlag):
	#----------------------------------------------------
	#定数の代入
	#----------------------------------------------------
	NS = 2000
	IS = 0.0005
	Xs0 = 35.0
	#----------------------------------------------------
	#ファイルの有無の確認
	#----------------------------------------------------
	if os.path.exists(filename):
		a = 1
	else:
		print(filename+"が存在しません")
		a = input()
		return
	#----------------------------------------------------
	#csvを開いて入力データを読み込む
	#----------------------------------------------------
	with open(filename,mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		E = []
		Amp  = 0.0
		Amp2 = 1.0
		if AmpFlag == 1: #歪増幅を考慮する場合は塑性率の換算係数Amp,Amp2を計算する
			Amp2= (1.0 / (LPL0 + (1.0-LPL0)*APAE))	#弾性域の増幅
			Amp = (APAE * ((1.0-LPL0)/LPL0))		#塑性域の増幅
			mucr = 1.0 / Amp2
		elif AmpFlag == 2: #歪増幅を考慮する場合は塑性率の換算係数Amp,Amp2を計算する
			Amp2= (1.0 / (LPL0 + (1.0-LPL0)*APAE))	#弾性域の増幅
			Amp = (1.0-LPL0)/LPL0					#塑性域の増幅
			mucr = 1.0 / Amp2
		for line in csvline:
			mu = float(line[1])
			if AmpFlag >= 1:
				if mu >= mucr:  mu = (mu*(1.0+Amp)+1.0-(1.0+Amp)/Amp2)*Ey #塑性域
				else:           mu = (mu*Amp2)*Ey          				  #弾性域
			E.append([float(line[0]),mu])
	#----------------------------------------------------
	#あらかじめ極値を分解したリストをつくる
	#----------------------------------------------------
	EE,IP,AS = FuncKokuchi(E)
	#----------------------------------------------------
	#歪み度範囲と破断繰返し回数のリスト作成
	#----------------------------------------------------
	S1 = IS * -0.5			#ある歪度範囲の中央値の計算に使う。
	S2 = (1.0+1.0) * -Ey	#全振幅で見たときは2倍(正+負)の降伏歪を取り除く。
	FR = []
	FR1 = []
	dum = S1 + S2
	Nf = []
	KP = 0.0
	MI = 0
	for i in range(NS):
		dum = dum + IS
		if MI == 0 and dum > 0.0:MI = i
		FR.append(dum)
		if FR[-1] < 0.0:FR[-1] = 0.0	#レインフロー法で統計した歪度数から弾性分を除いた塑性歪の集計値を出すために使う。
		KP = KP + IS
		FR1.append(round(KP,4))			#レインフロー法で統計する時の歪度範囲は弾性分を含んだ単純な形式
		if FR[i] > 0.0:
			Nf.append(FuncCoffinMansonNf(FR1[i],m2,C2*0.50))#マイナー則のための処理
		else:
			Nf.append(999999999999.0)						#マイナー則のための処理
	#----------------------------------------------------
	#フラグの初期化
	#----------------------------------------------------
	YieldFlag = 0
	FractureFlag = 0
	#----------------------------------------------------
	#損傷度時刻歴データファイルの出力準備
	#----------------------------------------------------
	if OutFlag == 1:
		file = open('Out_DamageHistory_'+filename,mode='w',newline="")
		out = csv.writer(file)
		out.writerow(['Time(s)','軸ひずみ(無次元)x増幅','実効の平均塑性ひずみ片振幅(%)','実効の累積塑性ひずみ(%)','絶対値最大軸ひずみ(%)','骨格比','疲労破壊が生じる累積塑性ひずみ(%)','判定(1:疲労破壊)','マイナー則で評価した時の損傷度D'])  
	#----------------------------------------------------
	#変数の初期化
	#----------------------------------------------------
	CumEp = 0.0				#実効の累積塑性ひずみ
	FinalAveAmp = 0.0 		#最終(または破断時)平均塑性ひずみ振幅
	FinalAs = 0.0 			#最終(または破断時)骨格比
	FinalFractureCumEp = 0.0#最終(または破断時)破断判定歪振幅
	#----------------------------------------------------
	#各ステップの損傷度評価
	#----------------------------------------------------
	for i in range(IP):
		if i <= 2: continue		#レインフロー法は3点以上のデータが必要
		#----------------------------------------------------
		#変数の初期化
		#----------------------------------------------------
		CumAmp = 0.0			#平均塑性歪を出すための見かけの累積塑性ひずみ
		AveAmp = 0.0			#実効の平均塑性ひずみ
		FratureCumEp = 0.0		#疲労破断判定用の累積塑性ひずみ
		Dminer = 0.0			#マイナー則評価用の損傷度D
		As = 0.0				#i時点の骨格比
		#----------------------------------------------------
		#平均塑性ひずみ片振幅の評価
		#----------------------------------------------------
		IFM = RainFlow(EE[0:i+1],i+1,NS,FR1)					#レインフロー法で頻度数を計数する
		Freq = sum(IFM[MI:]) 									#塑性域の頻度数を集計する
		for j in range(NS): CumAmp = CumAmp + FR[j]*IFM[j]		#頻度数に基づく見かけの累積塑性歪を集計する
		if Freq != 0.0: AveAmp = CumAmp / Freq * 100.0 * 0.5	#平均塑性ひずみ(片振幅)
		#----------------------------------------------------
		#累積塑性ひずみ振幅の評価
		#----------------------------------------------------
		if abs(EE[i][1]) >= Ey:YieldFlag = 1					#初期降伏の判定
		if i == 0:												#1ステップ目だけの特別な処理
			CumEp = (EE[i][1] - Ey) * 100.0
			if CumEp < 0.0:CumEp = 0.0
		elif YieldFlag == 1:									#初期降伏後はひずみ増分を累積塑性ひずみに累加する
			CumEp = CumEp + abs(EE[i][1] - EE[i-1][1])*100.0
		#----------------------------------------------------
		#骨格比の評価
		#----------------------------------------------------	
		if AveAmp > 0 and CumEp > 0.0:
			As = AS[i] - Ey
			if As < 0.0:As = 0.0
			As = AS[i] / CumEp #絶対値経験最大ひずみ / 累積塑性ひずみ振幅
		#----------------------------------------------------
		#疲労寿命の評価
		#----------------------------------------------------	
		FratureCumEp = As / Xs0 + ((1.0-As)*0.25)*((AveAmp**(1.0+m2))/C2)**(-1.0/m2)	#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ(分母を計算)
		if FratureCumEp > 0.0:FratureCumEp = 1.0 / FratureCumEp							#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ(逆数をとって完成)
		else: FratureCumEp = 99999999999.9												#分母が0.0以下は極大値をいれておく
		#----------------------------------------------------
		#疲労破壊の判定(破断が評価された場合はその時点の評価値を保存する)
		#----------------------------------------------------	
		if FratureCumEp <= CumEp:
			FractureFlag = 1
			FinalAveAmp = AveAmp
			FinalAs = As
			FinalFractureCumEp = FratureCumEp
		#----------------------------------------------------
		#マイナー則の損傷度評価(おまけ)
		#----------------------------------------------------
		for j in range(NS):
			if FR1[j] > 0.0: Dminer = Dminer + IFM[j] / Nf[j]	#弾性範囲も含めて淡々と頻度数と破断繰返し回数で判定
		#----------------------------------------------------
		#損傷度時刻歴を出力する場合の処理
		#----------------------------------------------------
		if OutFlag == 1:
			if FratureCumEp == 99999999999.9:FratureCumEp='N/A'
			out.writerow([EE[i][0],EE[i][1],AveAmp,CumEp,AS[i],As,FratureCumEp,FractureFlag,Dminer])
	if OutFlag:file.close()
	#----------------------------------------------------
	#出力用の処理
	#----------------------------------------------------
	if FractureFlag > 0:
		AveAmp       = FinalAveAmp
		FratureCumEp = FinalFractureCumEp
		As           = FinalAs
	elif FratureCumEp == 99999999999.9:FratureCumEp='N/A'
	#----------------------------------------------------
	#最終的な評価値をreturn
	#----------------------------------------------------
	return [AS[i],AveAmp,CumEp,FratureCumEp,As,FractureFlag,Dminer]

if os.path.exists('InputFileIndex.csv'):
	OutList = []
	with open('InputFileIndex.csv',mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		i = 0
		for i,line in enumerate(csvline):
			if i <= 0:continue      #FuncFractureBRB(filename,     AmpFlag,            Ey,          LPL0,          APAE,            m2,            C2,     OutFlag)
			OutList.append([line[0]]+FuncFractureBRB(line[0] ,int(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),int(line[7])))
			print(line[0]+"　評価終了")
	file = open('Out_DamageEvaluationBRB.csv',mode='w',newline="")
	out = csv.writer(file)
	out.writerow(['ファイル名','絶対値最大軸ひずみ(%)','最終(または破断時)の平均塑性ひずみ片振幅(%)','最終の累積塑性ひずみ(%)','最終(または破断時)の疲労破壊が生じる累積塑性ひずみ(%)','最終(または破断時)の骨格比','判定(1:疲労破壊)','マイナー則で評価した時の損傷度D'])
	for line in OutList:
		out.writerow(line)
	file.close()
else:
	print('InputFileIndex.csv'+"が存在しません")
	a=input()