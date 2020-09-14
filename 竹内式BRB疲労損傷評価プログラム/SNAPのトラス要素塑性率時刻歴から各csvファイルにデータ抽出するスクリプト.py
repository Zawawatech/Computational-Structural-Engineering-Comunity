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
#
#----------------------------------------------------
print("csvファイル名を記入せよ(拡張子込)")
Filename = input()
# Filename = "QDList.csv"
print("出力ファイルの頭の名前を入力せよ")
OutHead = input()
# OutHead = "Kobe"
with open(Filename,mode='r') as csvfile:
	csvline = csv.reader(csvfile)
	i = 0
	for i,line in enumerate(csvline):
		if   i > 2:
			k = 0
			for j in range(mem):
				out[j].writerow([line[k],line[k+1]])
				k = k + 2
		elif i == 0:continue
		elif i == 1:
			OutFile = []
			out = []
			mem = 0
			for dum in line:
				if dum != '':
					mem = mem+1
					print(OutHead+'_'+dum.strip()+'.csv')
					OF = open(OutHead+'_'+dum.strip()+'.csv',mode='w',newline="")
					OutFile.append(OF)
					out.append(csv.writer(OF))
		elif i == 2:continue
	for j in range(mem):OutFile[j].close()
