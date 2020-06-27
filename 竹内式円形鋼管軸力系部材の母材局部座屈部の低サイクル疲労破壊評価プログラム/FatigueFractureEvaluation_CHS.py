#-----------------------------------------------------
#　　繰返し変形による円形鋼管軸力部材の低サイクル疲労破壊
#									CODED by Y.T(2020)
#-----------------------------------------------------
#更新履歴
#2020/06/20		戯れに作る。
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
import math
#----------------------------------------------------
#Coffin-Manson則の破断繰り返し回数を計算する。
#引数：塑性ひずみ振幅,m2,C2
#戻り値：その振幅の破断繰り返し回数Nf
#----------------------------------------------------
def FractureEPRange(AveEpRange,C,m):
	FEPR = 2.0*C**(-1.0/m)*AveEpRange**((m+1)/m)
	return FEPR
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
def RainFlow(E,LY,NS,FR1,FR2):
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
				KP = 0.0
				for j in range(NS):
					if FM <= FR2[0]:	#1%にも満たない局部塑性歪は計数しないで吹き飛ばす。(2015以前の旧仕様)
						IP = IP - 2
						for l in range(AA,IP-2):
							Y[l] = Y[l+2]
						k=3
						break
					if FM <= FR1[j]:
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
			# if FM <= FR2[j]: continue
			if FM <= FR1[j]:
				IFM[j] = IFM[j] + 1.0
				break
	#---------------------------------------------------
	#Return
	#---------------------------------------------------
	return IFM
#----------------------------------------------------
#全体座屈発生歪の計算
#引数：D,t,lk,Sy
#戻り値：全体座屈発生歪(無次元)
#----------------------------------------------------
# def FuncEB(D,t,lk,Sy):
# 	d = D-2*t
# 	iy = (D**2.0+d**2.0)**0.5 / 4.0
# 	lamda = lk / iy
# 	Limitlamda = 3.141592653589*(205000.0/0.6/Sy)**0.5
# 	ro = (lamda / Limitlamda)**2.0
# 	# if lamda <= Limitlamda:Sb = ((1.0-0.4*ro)/(1.5+2.0/3.0*ro)*Sy) * 1.5	#鋼構造基準fcの1.5倍
# 	# if lamda >  Limitlamda:Sb = (0.277/ro*Sy) * 1.5						#鋼構造基準fcの1.5倍
# 	if lamda <= Limitlamda:Sb = ((1.0-0.4*ro)*Sy)      						#単純に安全率なし
# 	if lamda >  Limitlamda:Sb = (0.6/ro*Sy)      							#単純に安全率なし
# 	EB = (Sb/205000.0)
# 	return EB
#----------------------------------------------------
#局部座屈発生歪の計算
#引数：D,t,Sy
#戻り値：局部座屈発生歪(無次元)
#----------------------------------------------------
def FuncELB(D,t,Sy):
	Ey  = Sy / 205000.0
	ELB = 0.0683*(Ey**-0.39)*((D/t)**-1.39)
	return ELB
#----------------------------------------------------
#局部歪の計算
#引数：
#戻り値：局部歪(無次元)
#----------------------------------------------------
def LocalE(E,D,t,lk,Ey,EB,ELB,Time_B,Time_LB,Entm,CumEp,LERANGEPRE,Epre,time):
	def FuncBetaC(D,t):
		RDT = D/t
		if   RDT < 35.0:BetaC = 1.0
		elif RDT < 50.0:BetaC = (RDT-20.0)/15.0
		else           :BetaC = 2.0
		return BetaC
	RDT   = D/t
	BetaC = FuncBetaC(D,t)
	QLB   = math.acos(1.0-ELB)
	PAI   = 3.141592653589
	LP    = PAI * math.sqrt(D*t/6.0)
	#----------------------------------------------------
	#経験最大引張ひずみの更新
	#----------------------------------------------------	
	if Entm < E:Entm = E
	#----------------------------------------------------
	#等価軸ひずみ振幅と塑性回転角QHの計算
	#----------------------------------------------------
	EQRANGE = Entm - E
	QH      = math.acos(1.0-EQRANGE)
	#----------------------------------------------------
	#座屈状態に基づいた局部塑性ひずみ振幅の評価
	#----------------------------------------------------
	if Time_LB != 0.0 and time>Time_LB:#強制的に局部座屈とする場合。
		dQH    = QH-QLB
		FAI    = math.cos(dQH)-D*math.sin(dQH)/LP
		if FAI <= -1.0: FAI = PAI
		else:           FAI = math.acos(FAI)
		LE     = -(3.0*(6.0)**(0.5)*FAI*BetaC/(2.0*PAI*(RDT)**(0.5)) + (D*QLB)/(lk*(1-PAI*0.25)))
		LERANGE= Entm + LE
		AlphaC = LE / EQRANGE * -1.0
		flag   = 2
	elif Time_B != 0.0 and time > Time_B:#強制的に全体座屈させる場合
		LE     = E
		LERANGE= EQRANGE
		AlphaC = 1.0
		flag   = 0
	elif EQRANGE <= EB or time < Time_B :#全体座屈前
		LE     = E
		LERANGE= EQRANGE
		AlphaC = 1.0
		flag   = 0
	elif EQRANGE < ELB or time < Time_LB:#全体座屈後、局部座屈前
		LE     = -(QH*D) / (lk*(1.0-PAI*0.25))
		AlphaC = LE / EQRANGE * -1.0
		LERANGE= Entm + LE
		flag   = 1
	else:           #局部座屈後
		dQH    = QH-QLB
		FAI    = math.cos(dQH)-D*math.sin(dQH)/LP
		if FAI <= -1.0: FAI = PAI
		else:           FAI = math.acos(FAI)
		LE     = -(3.0*(6.0)**(0.5)*FAI*BetaC/(2.0*PAI*(RDT)**(0.5)) + (D*QLB)/(lk*(1-PAI*0.25)))
		LERANGE= Entm + LE
		AlphaC = LE / EQRANGE * -1.0
		flag   = 2
	#----------------------------------------------------
	#累積塑性歪の評価
	#----------------------------------------------------
	if time == 0:
		DifEp = 0.0
		CumEp = (LERANGE-E)*100.0
	else:
		DifEp = math.fabs(LERANGE-LERANGEPRE) + (Epre-E)
		CumEp = CumEp + DifEp*100.0
	LERANGEPRE = LERANGE
	Epre = E
	#----------------------------------------------------
	#リターン
	#----------------------------------------------------
	return EQRANGE,LE,LERANGE,LERANGEPRE,Epre,CumEp,DifEp,AlphaC,Entm,flag
#----------------------------------------------------
#破断評価のターミナルルーチン
#引数：軸ひずみ時刻歴ファイル名,D,t,座屈長さlk,降伏応力度Sy,座屈応力度SbC,m,全体座屈発生時間,局部座屈発生時間,全体座屈発生歪,局部座屈発生歪
#戻り値：
#----------------------------------------------------
def FuncFractureCHS(filename,D,t,lk,Sy,Sb,C,m,Time_B,Time_LB,E_B,E_LB,OutFlag):
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
			E.append([float(line[0]),float(line[1])])
	#----------------------------------------------------
	#部材諸元の計算
	#----------------------------------------------------
	if E_B  > 0.0:EB = E_B
	else:         EB = Sb/205000.0
	if E_LB > 0.0:ELB= E_LB
	else:         ELB= FuncELB(D,t,Sy)
	Ey = Sy/205000.0
	#----------------------------------------------------
	#局部歪の評価
	#----------------------------------------------------
	EP  = []
	CEP = []
	dEP = []
	Ac  = []
	Et  = []
	dE  = []
	flag = []
	dEp = []
	dLE = []
	Entm = 0.0
	Epre  = 0.0
	CumEp = 0.0
	DifEp = 0.0
	LERANGEPRE = 0.0
	for line in E:
		EQRANGE,LE,LERANGE,LERANGEPRE,Epre,CumEp,DifEp,AlphaC,Entm,dflag = LocalE(line[1],D,t,lk,Ey,EB,ELB,Time_B,Time_LB,Entm,CumEp,LERANGEPRE,Epre,line[0])
		EP.append([line[0],LERANGE])
		dE.append(EQRANGE)
		CEP.append(CumEp)
		dEp.append(DifEp)
		Ac.append(AlphaC)
		flag.append(dflag)
		dLE.append(LE)
	#----------------------------------------------------
	#極値分解
	#----------------------------------------------------
	EE,IP = FuncKokuchi(EP)
	#----------------------------------------------------
	#歪み度範囲の作成
	#----------------------------------------------------
	IS = 0.02
	NS = 50
	S1 = IS * -0.5			#ある歪度範囲の中央値の計算に使う。
	S2 = EB+Ey          	#全振幅で見たときは2倍(正+負)の降伏歪を取り除く。
	FR = []
	FR1 = []
	FR2 = []
	dum = S1 + S2
	KP = 0.0
	for i in range(NS):
		dum = dum + IS
		FR.append(dum)
		if FR[-1] < 0.0:FR[-1] = 0.0	#レインフロー法で統計した歪度数から弾性分を除いた塑性歪の集計値を出すために使う。
		KP = KP + IS
		FR1.append(round(KP,4))			#レインフロー法で統計する時の歪度範囲は弾性分を含んだ単純な形式
		FR2.append(round(KP*0.5,4))
	#----------------------------------------------------
	#フラグの初期化
	#----------------------------------------------------
	FractureFlag = 0
	#----------------------------------------------------
	#中間出力ファイルの準備
	#----------------------------------------------------
	if OutFlag == 1:
		file = open('Out_DamageHistory_'+filename,mode='w',newline="")
		out = csv.writer(file)
		out.writerow(['Time(s)','等価軸ひずみ(%)','等価軸ひずみ振幅(%)','歪振幅拡大係数Ac','局部歪(%)','局部歪振幅(%)','平均塑性歪振幅(%)','累積塑性歪(%)','対応する亀裂発生の累積塑性ひずみ(%)','弾性(0),全体座屈(1),局部座屈(2),破断(3)'])
	#----------------------------------------------------
	#各ステップの損傷度評価
	#----------------------------------------------------
	IIP = -1
	FractureFlag = 0
	for i in range(len(E)):
		print(filename,E[i][0])
		AveEpRange = 0.0			#実効の平均塑性ひずみ
		CumRange = 0.0				#平均塑性歪振幅を出すための見かけの累積塑性ひずみ
		FratureCumEp = 999999999.99	#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ
		dflag = 0
		if IIP+1 <= len(EE)-1:
			if EP[i][0] == EE[IIP+1][0]:IIP = IIP + 1
		if IIP > 2:
			IFM = RainFlow(EE[0:IIP+1],IIP+1,NS,FR1,FR2)
			Freq = sum(IFM)
			for j in range(NS):CumRange = CumRange + FR[j]*IFM[j]	#頻度数に基づく累積歪
			if Freq != 0.0:
				AveEpRange   = CumRange / Freq * 100.0 	#平均塑性ひずみ(全振幅)
				FratureCumEp = FractureEPRange(AveEpRange,C,m)			#現時点の平均塑性ひずみに対応した破断判定の累積塑性ひずみ
		if FratureCumEp <= CEP[i]: FractureFlag = 3					#破断判定
		if OutFlag == 1:
			dflag = flag[i]
			if FractureFlag == 3: dflag = FractureFlag
			out.writerow([E[i][0],E[i][1]*100.0,dE[i]*100.0,Ac[i],dLE[i]*100.0,EP[i][1]*100.0,AveEpRange,CEP[i],FratureCumEp,dflag])
	if OutFlag:file.close()
	#----------------------------------------------------
	#return
	#----------------------------------------------------
	return [AveEpRange,CEP[-1],FratureCumEp,FractureFlag]

if os.path.exists('InputFileIndex.csv'):
	OutList = []
	with open('InputFileIndex.csv',mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		i = 0
		for i,line in enumerate(csvline):
			if i <= 0:continue
			OutList.append([line[0]]+FuncFractureCHS(line[0],
				float(line[1]),float(line[2]),float(line[3]),
				float(line[4]),float(line[5]),float(line[6]),float(line[7]),
				float(line[8]),float(line[9]),float(line[10]),
				float(line[11]),int(line[12])))
			print(line[0]+"　評価終了")
	file = open('Out_DamageEvaluationCHS.csv',mode='w',newline="")
	out = csv.writer(file)
	out.writerow(['ファイル名','最終の平均塑性ひずみ振幅(%)','最終の累積塑性ひずみ(%)','疲労破壊が生じる累積塑性ひずみ(%)','判定(1:疲労破壊)'])
	for line in OutList:
		out.writerow(line)
	file.close()
else:
	print('InputFileIndex.csv'+"が存在しません")
	a=input()