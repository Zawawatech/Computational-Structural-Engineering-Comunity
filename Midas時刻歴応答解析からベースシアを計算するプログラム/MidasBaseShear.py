#-----------------------------------------------------
#　　繰返し変形による座屈拘束ブレース(BRB)疲労損傷の検証方法
#									CODED by Y.T(2019)
#-----------------------------------------------------
#更新履歴
#2019/08/29		学生の要望に答えて作る。
#-----------------------------------------------------
#-----------------------------------------------------
#必要モジュールのimport
#-----------------------------------------------------
import sys         	#OSのコマンドや環境変数を使う
import csv         	#csvファイルの処理
import os.path     	#OSのコマンド
import math

#-----------------------------------------------------
#部材せん断力を全体座標系に変更した上で総和の時刻歴を計算する。
#-----------------------------------------------------
#Infile:入力ファイル名
#DT:時間刻み
#Angle:コードアングル
#OutFile:出力ファイル名
#-----------------------------------------------------
def SHEAR(Infile,DT,Angle,OutFile):
	#-----------------------------------------------------
	#ファイルを開く
	#-----------------------------------------------------
	f = open(Infile,mode='r')
	dum = 0
	flag = 0
	Time = 0.0
	Inimem = 0
	RowNum = 0
	TotalForce = []
	#-----------------------------------------------------
	#角度をラジアンに変換する。
	#-----------------------------------------------------
	Ag = math.pi * float(Angle) / 180.0
	#-----------------------------------------------------
	#中身を読み込んで総和する。
	#-----------------------------------------------------
	while dum == 0:
		line = f.readline()
		if "要素" in line: print(line)
		if not line:
			dum = 1
			continue
		if " 選択した梁要素の最大/最小 断面力       単位系 ..: mm, kN, sec" in line:
			dum = 1
			continue
		if   flag == 0 and "  TIME        AXIAL       SHEAR-y      SHEAR-z      TORSION      MOMENT-y     MOMENT-z " in line:
			flag = 1
			continue
		elif flag == 1 and "------- -  -----------  -----------  -----------  -----------  -----------  ----------- " in line:
			flag = 2
			RowNum = 0
			continue
		elif flag == 2 and "J" in line:
			Fx = float(line[ 9:22])
			Fy = float(line[22:35])
			Fz = float(line[35:48])
			Force = [Time,Fz*math.cos(Ag) + Fy*math.sin(Ag),Fz*math.sin(Ag) - Fy*math.cos(Ag),Fx]
			if Inimem == 0:
				TotalForce.append(Force)	#初めての部材はそのまま足し込む
			else:							#2番目以降の部材はTotalForceに足し込んでいく
				for i in range(1,4):
					TotalForce[RowNum][i] = TotalForce[RowNum][i] + Force[i]
			Time = Time + DT
			RowNum = RowNum + 1
		elif flag == 2  and "------- -  -----------  -----------  -----------  -----------  -----------  ----------- " in line:
			flag = 0
			Time = 0.0
			Inimem = 1
			continue
	f.close()
	#-----------------------------------------------------
	#書き出す
	#-----------------------------------------------------
	Outt = open(OutFile,mode='w',newline="")
	out = csv.writer(Outt)
	out.writerow(['time(s)','Base-Fx(kN)','Base-Fy(kN)','Base-Fz(kN)'])
	for line in TotalForce:
		out.writerow(line)
	Outt.close()

print("Midasから出力された柱の断面力時刻歴(部材座標系)のテキストファイル名(拡張子をふくむ)を入力して下さい。")
filename = input()
print("出力の時間刻み(s)を入力して下さい。")
dt = float(input())
print("全体座標Z軸に対して柱のコードアングル(°)を入力して下さい。")
Angle = float(input())
SHEAR(filename,dt,Angle,'Out.csv')
