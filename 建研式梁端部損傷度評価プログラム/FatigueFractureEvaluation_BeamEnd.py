#-----------------------------------------------------
#　　　　　　繰返し変形による梁端部破断の検証方法
#									CODED by Y.T(2020)
#-----------------------------------------------------
#更新履歴
#2020/05/28		戯れに作る。
#-----------------------------------------------------
#----------------------------------------------------
#ModuleのImport
#----------------------------------------------------
import csv
import os
#----------------------------------------------------
#梁端の降伏回転角を計算するサブルーチン
#引数：ファイル名,θP,MP,時刻歴データの出力フラグ
#戻り値：損傷度、最終塑性率、最終累積塑性変形倍率、θp、Mp
#----------------------------------------------------
def FuncIta(filename,C,Beta,SitaP,MP,OutFlag):
	#----------------------------------------------------
	#ファイルの有無の確認
	#----------------------------------------------------	
	if os.path.exists(filename):
		a = 1
	else:
		print(filename+"が存在しません")
		a = input()
		return 0.0,0.0,0.0,0.0,0.0
	with open(filename,mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		SM = []
		for line in csvline:
			SM.append([float(line[0]),float(line[1])])
	#----------------------------------------------------
	#初期化
	#----------------------------------------------------
	Wp = 0.0
	SitaMax = 0.0
	Wplist = [0.0]
	SitaMaxlist = [0.0]
	Italist = [0.0]
	Dlist = [0.0]
	#----------------------------------------------------
	#荷重変形関係から諸量を計算する。
	#----------------------------------------------------
	for num,line in enumerate(SM):
		if num==0:
			line2 = line
			continue
		dS = line[0]-line2[0]		#回転角増分
		dM = (line[1]-line2[1])  	#前後のモーメントの平均値
		if MP == 0.0 or SitaP == 0.0:						#MpとSitaPの計算
			if dS != 0.0 and abs(dM/dS) <= 0.01:
				if SitaP == 0:SitaP = abs(line2[0])
				if MP == 0:MP = abs(line2[1])
		if SitaMax < abs(line[0]):SitaMax = abs(line[0])	#最大回転角の更新
		Wp = Wp + (line[1]+line2[1])*dS*0.5					#累積吸収エネルギーの更新(区分求積法)
		line2 = line
		if OutFlag==1:
			Wplist.append(Wp)
			SitaMaxlist.append(SitaMax)
	#----------------------------------------------------
	#最終損傷度の計算
	#----------------------------------------------------
	Ita2 = Wp / (MP*SitaP)
	Mu2 = SitaMax / SitaP
	Damage2 = Ita2 / (4.0*(Mu2 - 1.0)) * (Mu2 / C) **(1.0 / Beta)
	#----------------------------------------------------
	#必要に応じて損傷度時刻歴を計算する。
	#----------------------------------------------------
	if OutFlag == 1:
		file = open('Out_DamageHistory_'+filename,mode='w',newline="")
		out = csv.writer(file)
		out.writerow(['塑性率','累積塑性変形倍率','損傷度D'])  
		for i in range(num):
			Wp = Wplist[i]
			SitaMax = SitaMaxlist[i]
			Ita = Wp / (MP*SitaP)
			Mu = SitaMax / SitaP
			Damage = 0.0
			if Mu > 1.0:Damage = Ita / (4.0*(Mu - 1.0)) * (Mu / C) **(1.0 / Beta)
			out.writerow([Mu,Ita,Damage])
		file.close()
	#----------------------------------------------------
	#Return
	#----------------------------------------------------		
	return [Mu2,Ita2,Damage2,SitaP,MP]

if os.path.exists('InputFileIndex.csv'):
	OutList = []
	with open('InputFileIndex.csv',mode='r') as csvfile:
		csvline = csv.reader(csvfile)
		i = 0
		for i,line in enumerate(csvline):
			if i <= 0:continue
			OutList.append([line[0]]+FuncIta(line[0],float(line[1]),float(line[2]),float(line[3]),float(line[4]),int(line[5])))
	file = open('Out_DamageEvaluation.csv',mode='w',newline="")
	out = csv.writer(file)
	out.writerow(['ファイル名','塑性率','累積塑性変形倍率','損傷度D','θp','Mp'])
	for line in OutList:
		out.writerow(line)
	file.close()
else:
	print('InputFileIndex.csv'+"が存在しません")
	a=input()