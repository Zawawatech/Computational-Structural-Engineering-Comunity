#-----------------------------------------------------
#　　繰返し変形による鋼材の低サイクル疲労損傷の評価プログラム
#									CODED by Y.T(2020)
#-----------------------------------------------------
#更新履歴
#2020/06/07		戯れに作る。※塑性ひずみの計算から弾性ひずみは消去しています。
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
#破断繰り返し回数を計算する。
#引数：振幅,C,m
#戻り値：その振幅の破断繰り返し回数Nf
#----------------------------------------------------
def FuncNf(dEp,C,m):
	Nf = C*(dEp) ** (m)
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
#引数：時刻歴ファイル名,マイナー則の単位振幅サイズ,マイナー則の度数,評価タイプ1,評価タイプ2,塑性変形,C,m2,損傷度時刻歴の出力フラグ
#戻り値：損傷度D,破断判定
#----------------------------------------------------
def FuncFractureBRB(filename,IS,NS,EveType1,EveType2,Ey,C,m,OutFlag):
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
	#極値分解
	#----------------------------------------------------
	EE,IP = FuncKokuchi(E)
	#----------------------------------------------------
	#歪み度範囲の作成
	#----------------------------------------------------
	S1 = 0.0
	S2 = 0.0
	if EveType1 == 1:IS = IS*2.0			#片振幅評価のときは一旦2倍にする。
	S1 = IS * -0.5							#ある歪度範囲の中央値の計算に使う。
	if EveType2 == 1:S2 = (1.0+1.0) * -Ey	#全振幅で見たときは2倍(正+負)の降伏歪を取り除く。
	FR = []
	FR1 = []
	FR2 = []
	dum = S1 + S2
	Nf = []
	KP = 0.0
	for i in range(NS):
		dum = dum + IS
		FR.append(dum)
		if FR[-1] < 0.0:FR[-1] = 0.0	#レインフロー法で統計した歪度数から弾性分を除いた塑性歪の集計値を出すために使う。
		KP = KP + IS
		FR1.append(round(KP,4))			#レインフロー法で統計する時の歪度範囲は弾性分を含んだ単純な形式
		FR2.append(round(KP*0.5,4))
		if FR[i] > 0.0:
			if EveType1 == 0:Nf.append(FuncNf(FR[i],C,m))		#全振幅評価時
			else:            Nf.append(FuncNf(FR[i]*0.5,C,m))	#片振幅評価時
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
		out.writerow(['Time(s)','変形量(入力値)','実効の累積損傷度D','判定(1:疲労破壊)']+FR2)  
	#----------------------------------------------------
	#各ステップの損傷度評価
	#----------------------------------------------------
	for i in range(IP):
		if i <= 2: continue
		IFM = RainFlow(EE[0:i+1],i+1,NS,FR1)
		Freq = sum(IFM)
		Dminer = 0.0	#マイナー則評価用の損傷度D
		for j in range(NS):
			if FR[j] > 0.0: Dminer = Dminer + IFM[j] / Nf[j]
			if Dminer > 1.0: FractureFlag = 1.0
		if OutFlag == 1:
			out.writerow([EE[i][0],EE[i][1],Dminer,FractureFlag]+[int(j) for j in IFM])
	if OutFlag:file.close()
	#----------------------------------------------------
	#return
	#----------------------------------------------------
	return [Dminer,FractureFlag]

if os.path.exists('InputFileIndex.csv'):
	OutList = []
	with open('InputFileIndex.csv',mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		i = 0
		for i,line in enumerate(csvline):
			if i <= 0:continue
			OutList.append([line[0]]+FuncFractureBRB(line[0],float(line[1]),int(line[2]),int(line[3]),int(line[4]),float(line[5]),float(line[6]),float(line[7]),int(line[8])))
			print(line[0]+"　評価終了")
	file = open('Out_DamageEvaluation.csv',mode='w',newline="")
	out = csv.writer(file)
	out.writerow(['ファイル名','最終の累積疲労損傷度D','判定(1:疲労破壊)'])
	for line in OutList:
		out.writerow(line)
	file.close()
else:
	print('InputFileIndex.csv'+"が存在しません")
	a=input()