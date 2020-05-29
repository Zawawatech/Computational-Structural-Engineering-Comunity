!********************************************************
!   Subroutine  :   Driver_PlasticityModel_Uniaxial_YamashitaPinBearing
!   Purpose     :   To compute normalized material nonlinearity 
!                  (i.e. Inner force, tangent modulas)
!                   無次元量として復元力特性を計算する。
!   Reference   :   Watanabe S. and Yamashita T.,
!                   Earthquake response analysis of steel roof gymnasiums considering nonlinear
!                   restoring force characteristics of lower structures and roof bearings,
!                   Trans. AIJ ("Kibyo-shi 黄表紙"), AIJ, Vol.85, No.768,pp.209-218, 2020.2.
!   Remarks     :   This subroutine is based on an updated backbone curve tecnique by Yuki TERAZAWA (2020)
!                   このルーチンは骨格曲線を逐次更新するテクニックを用いて書き下ろしました。
!               :   Rule4のピンチング前のスリップを表現するために特別な折れ点13,16を加えています。
!********************************************************    
subroutine Driver_PlasticityModel_Uniaxial_YamashitaPinBearing(&
            istep,iter,istage,&
            LD,LS,&
            D,S,T,&
            DD,&
            P1,P2,P3,P13,P4,P5,P6,P16,&
            Q1,Q2,Q3,Q13,Q4,Q5,Q6,Q16,&
            T1,T2,T3,T6,T7,&
            IP1,IP2,IP3,IP4,IP5,IP6,&
            IQ1,IQ2,IQ3,IQ4,IQ5,IQ6,&
            IT1,IT2,IT3,&
            nn,mu,C1,C2,C3,C4,C5,&
            Pmax,Pmin,Qmax,Qmin,Gyp,Gyn)
!--------------------------------------------------------
!implicite none
!--------------------------------------------------------        
    implicit none
!--------------------------------------------------------
!Define input & output arguments
!入力変数を宣言する
!--------------------------------------------------------
![in]
!   [integer]
!       istep   : Step number (今のステップ番号)    
!       iter    : Iteration number in this step (今の収斂回数番号)
!   [double precision]
!       DD      : The normalized incremental deformation of n-th iteration in this step (今回の無次元化増分変形)
!       IP1     : (Positive deformation) Break point 1 in the initial backbone curve (初期の骨格曲線折れ点1,正側の滑り出し,変形)
!       IP2     : (Positive deformation) Break point 2 in the initial backbone curve (初期の骨格曲線折れ点2,正側のBPLとアンカーボルト接触,変形)
!       IP3     : (Positive deformation) Break point 3 in the initial backbone curve (初期の骨格曲線折れ点3,正側のアンカーボルト曲げ降伏,変形)
!       IP4     : (Negative deformation) Break point 4 in the initial backbone curve (初期の骨格曲線折れ点4,負側の滑り出し,変形)
!       IP5     : (Negative deformation) Break point 5 in the initial backbone curve (初期の骨格曲線折れ点5,負側のBPLとアンカーボルト接触,変形)
!       IP6     : (Negative deformation) Break point 6 in the initial backbone curve (初期の骨格曲線折れ点6,負側のアンカーボルト曲げ降伏,変形)    
!       IQ1     : (Positive force) Break point 1 in the initial backbone curve (初期の骨格曲線折れ点1,正側の滑り出し,復元力)
!       IQ2     : (Positive force) Break point 2 in the initial backbone curve (初期の骨格曲線折れ点2,正側のBPLとアンカーボルト接触,復元力)
!       IQ3     : (Positive force) Break point 3 in the initial backbone curve (初期の骨格曲線折れ点3,正側のアンカーボルト曲げ降伏,復元力)
!       IQ4     : (Negative force) Break point 4 in the initial backbone curve (初期の骨格曲線折れ点4,負側の滑り出し,復元力)
!       IQ5     : (Negative force) Break point 5 in the initial backbone curve (初期の骨格曲線折れ点5,負側のBPLとアンカーボルト接触,復元力)
!       IQ6     : (Negative force) Break point 6 in the initial backbone curve (初期の骨格曲線折れ点6,負側のアンカーボルト曲げ降伏,復元力)        
!       IT1     : (Tangent modulas) Slope 1 in the initial backbone curve (初期の骨格曲線の勾配1,山下モデルのStage 1 )
!       IT2     : (Tangent modulas) Slope 2 in the initial backbone curve (初期の骨格曲線の勾配2,山下モデルのStage 2 )
!       IT3     : (Tangent modulas) Slope 3 in the initial backbone curve (初期の骨格曲線の勾配3,山下モデルのStage 3 )
!       nn      : the number of anchorbolt. (アンカーボルト本数)
!       mu      : friction coeficient (摩擦係数)
!       C1      : C1 = {Ny*Ny*(h+lm)/Mp}/Qf = N/Qf
!       C2      : C2 = {Mp/(h+lm)}/Qf
!       C3      : C3 = Pv/Qf
!       C4      : C4 = (h+lm)/δf
!       C5      : C5 = Ny/Qf
![inout]
!   [integer]    
!       istage  : ((n-1)-th iteration) The last stage  (前回のステージ番号)    
!   [double precision]
!       LD      : ((n-1)-th iteration) The last normalized deformation    (前回の無次元化変形)
!       LS      : ((n-1)-th iteration) The last normalized stress         (前回の無次元化荷重)
!       D       : (n-th iteration)     The current normalized deformation (今回の無次元化変形)
!       S       : (n-th iteration)     The current normalized stress      (今回の無次元化荷重)
!       T       : (n-th iteration)     The current tangent modulas        (今回の接線剛性比)
!       P1      : (Positive deformation) Break point 1 in the curent backbone curve (今のステップの骨格曲線折れ点1,正側の滑り出し,変形)
!       P2      : (Positive deformation) Break point 2 in the curent backbone curve (今のステップの骨格曲線折れ点2,正側のBPLとアンカーボルト接触,変形)
!       P3      : (Positive deformation) Break point 3 in the curent backbone curve (今のステップの骨格曲線折れ点3,正側のアンカーボルト曲げ降伏,変形)
!       P13     : (Positive deformation) Break point13 in the curent backbone curve (今のステップの骨格曲線折れ点13,再載荷時の正側のアンカーボルトのピンチング開始点,変形)
!       P4      : (Negative deformation) Break point 4 in the curent backbone curve (今のステップの骨格曲線折れ点4,負側の滑り出し,変形)
!       P5      : (Negative deformation) Break point 5 in the curent backbone curve (今のステップの骨格曲線折れ点5,負側のBPLとアンカーボルト接触,変形)
!       P6      : (Negative deformation) Break point 6 in the curent backbone curve (今のステップの骨格曲線折れ点6,負側のアンカーボルト曲げ降伏,変形)
!       P16     : (Negative deformation) Break point16 in the curent backbone curve (今のステップの骨格曲線折れ点16,再載荷時の正側のアンカーボルトのピンチング開始点,変形)
!       Q1      : (Positive force) Break point 1 in the curent backbone curve (今のステップの骨格曲線折れ点1,正側の滑り出し,復元力)
!       Q2      : (Positive force) Break point 2 in the curent backbone curve (今のステップの骨格曲線折れ点2,正側のBPLとアンカーボルト接触,復元力)
!       Q3      : (Positive force) Break point 3 in the curent backbone curve (今のステップの骨格曲線折れ点3,正側のアンカーボルト曲げ降伏,復元力)
!       Q13     : (Positive force) Break point13 in the curent backbone curve (今のステップの骨格曲線折れ点13,再載荷時の正側のアンカーボルトのピンチング開始点,復元力)
!       Q4      : (Negative force) Break point 4 in the curent backbone curve (今のステップの骨格曲線折れ点4,負側の滑り出し,復元力)
!       Q5      : (Negative force) Break point 5 in the curent backbone curve (今のステップの骨格曲線折れ点5,負側のBPLとアンカーボルト接触,復元力)
!       Q6      : (Negative force) Break point 6 in the curent backbone curve (今のステップの骨格曲線折れ点6,負側のアンカーボルト曲げ降伏,復元力)
!       Q16     : (Negative force) Break point16 in the curent backbone curve (今のステップの骨格曲線折れ点16,再載荷時の正側のアンカーボルトのピンチング開始点,復元力)   
!       T1      : (Tangent modulas) Slope 1 in the current backbone curve (今のステップの骨格曲線の勾配1,山下モデルのStage 1,正負,滑る前)
!       T2      : (Tangent modulas) Slope 2 in the current backbone curve (今のステップの骨格曲線の勾配2,山下モデルのStage 2,正側,滑り中)
!       T3      : (Tangent modulas) Slope 3 in the current backbone curve (今のステップの骨格曲線の勾配3,山下モデルのStage 3,正側,アンカー曲げ変形中)
!       T6      : (Tangent modulas) Slope 6 in the current backbone curve (今のステップの骨格曲線の勾配6,山下モデルのStage 2,負側,滑り中)
!       T7      : (Tangent modulas) Slope 7 in the current backbone curve (今のステップの骨格曲線の勾配7,山下モデルのStage 3,負側,アンカー曲げ変形中)
!       Pmax    : Maximum value of experienced deformation (経験最大変形)
!       Pmin    : Minimum value of experienced deformation (経験最小変形)    
!       Qmax    : Maximum value of experienced force (経験最大復元力)
!       Qmin    : Minimum value of experienced force (経験最小復元力)
!       Gyp     : ピンチングが始まる変形点(正側)
!       Gyn     : ピンチングが始まる変形点(負側)
!--------------------------------------------------------
    integer,            intent(   in)   :: iter,istep
    double precision,   intent(   in)   :: DD
    double precision,   intent(   in)   :: IP1,IP2,IP3,IP4,IP5,IP6
    double precision,   intent(   in)   :: IQ1,IQ2,IQ3,IQ4,IQ5,IQ6
    double precision,   intent(   in)   :: IT1,IT2,IT3
    double precision,   intent(   in)   :: nn,mu,C1,C2,C3,C4,C5
    integer,            intent(inout)   :: istage
    double precision,   intent(inout)   :: LD,LS
    double precision,   intent(inout)   :: D,S,T
    double precision,   intent(inout)   :: P1,P2,P3,P13,P4,P5,P6,P16
    double precision,   intent(inout)   :: Q1,Q2,Q3,Q13,Q4,Q5,Q6,Q16
    double precision,   intent(inout)   :: T1,T2,T3,T6,T7
    double precision,   intent(inout)   :: Pmin,Pmax,Qmin,Qmax,Gyp,Gyn
!--------------------------------------------------------
!Define local variables
!ローカル変数の宣言    
!--------------------------------------------------------
![integer]
!   stage   : Stage number in the backbone curve (骨格曲線内のステージ区分)
!--------------------------------------------------------    
    integer             :: stage
!--------------------------------------------------------
!Update backbone curve if iter = 0 
!収斂0回目なら骨格曲線を更新する。
!--------------------------------------------------------
    if(iter == 0) then
        if (istep == 1) then    !解析開始時点の初期値代入
            P1 = IP1
            P2 = IP2
            P3 = IP3
            P13 = IP3
            P4 = IP4
            P5 = IP5
            P6 = IP6
            P16 = IP6
            Q1 = IQ1
            Q2 = IQ2
            Q3 = IQ3
            Q13 = IQ3
            Q4 = IQ4
            Q5 = IQ5
            Q6 = IQ6
            Q16 = IQ6
            T1 = IT1
            T2 = IT2
            T3 = IT3
            T6 = IT2
            T7 = IT3
            Pmin = 0.0d0
            Pmax = 0.0d0
            Qmin = 0.0d0
            Qmax = 0.0d0
            Gyn = IP6
            Gyp = IP3
        else
            select case(istage) !前回のステップの最終収斂時点の荷重変形関係の座標から今回のステップの骨格曲線を作る。
                case(2) !Stage 2 (荷重が正側のサイド)
                    !----正側の骨格曲線は更新なし----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P2
                    Q2 = Q2
                    T2 = T2
                    P3 = P3
                    Q3 = Q3
                    T3 = T3
                    P13 = P13
                    Q13 = Q13
                    Gyp = Gyp
                    !----負側の骨格曲線は更新あり----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    if     (LD <  IP2 .and. LD > IP5)                   then !アンカーボルトとBPLが接触していない範囲, Rule 1 or Rule 3
                        Q5 = IQ5
                        P5 = IP5
                        T6 = (Q4-Q5) / (P4-P5)
                        if     (Pmin >  IP6) then !負側の曲げ降伏非経験, Rule 1
                            P6 = IP6
                            Q6 = IQ6
                            T7 = IT3
                            P16 = P6
                            Q16 = Q6
                        else if(Pmin <= IP6) then !負側の曲げ降伏　経験, Rule 3
                            T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                            Q6 = IQ6
                            P6 = P5 + (Q6-Q5)/T7
                            P16 = Pmin                    !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                            Q16 = Q6 + IT2 * (P16-P6) !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                        end if
                    else if(LD <= IP5 .and. Pmin >  IP6) then !アンカーボルトとBPLが接触している範囲, 曲げ降伏非経験, Rule 2
                        P5 = Pmin
                        Q5 = Qmin
                        T6 = (Q4 - Q5) / (P4 - P5)
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(LD <= IP5 .and. Pmin <= IP6) then !アンカーボルトとBPLが接触している範囲, 曲げ降伏　経験, Rule 4
                        Q5 = Q4 !負側のステージ2は一時的に消失する
                        P5 = P4 !負側のステージ2は一時的に消失する
                        T6 = T1 !負側のステージ2は一時的に消失する
                        P6 = Pmin
                        Q6 = IQ6
                        T7 = (Q5-Q6) / (P5-P6)
                        P16 = P6
                        Q16 = Q6
                    end if
                    Gyn = P16
                case(3) !Stage 3 (荷重が正側のサイド)
                    !----正側の骨格曲線は更新なし----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !正側のStage2は一時的に消失する
                    Q2 = Q1 !正側のStage2は一時的に消失する
                    T2 = T1 !正側のStage2は一時的に消失する
                    P3 = P3
                    Q3 = Q3
                    T3 = T3
                    P13 = P13
                    Q13 = Q13
                    Gyp = Gyp
                    !----負側の骨格曲線は更新あり----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    Q5 = IQ5
                    P5 = IP5
                    T6 = (Q4-Q5)/(P4-P5)
                    if     (Pmin >  IP6) then !負側の曲げ降伏非経験, Rule 1
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(Pmin <= IP6) then !負側の曲げ降伏　経験, Rule 3
                        T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                        Q6 = IQ6
                        P6 = P5 + (Q6-Q5)/T7
                        P16 = Pmin                    !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                        Q16 = Q6 + IT2 * (P16-P6) !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                    end if
                    Gyn = P16
                case(4) !Stage 4 (荷重が正側のサイド)
                    !----正側の骨格曲線は更新なし----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !正側のStage2は一時的に消失する
                    Q2 = Q1 !正側のStage2は一時的に消失する
                    T2 = T1 !正側のStage2は一時的に消失する
                    P3 = P1 !正側のStage3は一時的に消失する
                    Q3 = Q1 !正側のStage3は一時的に消失する
                    T3 = T1 !正側のStage3は一時的に消失する
                    P13 = P3
                    Q13 = Q3
                    Gyp = IP3 !正側のStage4の変形開始点は初期曲げ降伏点と同じ
                    !----負側の骨格曲線は更新あり----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    Q5 = IQ5
                    P5 = IP5
                    T6 = (Q4-Q5)/(P4-P5)
                    if     (Pmin >  IP6) then !負側の曲げ降伏非経験, Rule 1
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(Pmin <= IP6) then !負側の曲げ降伏　経験, Rule 3
                        T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                        Q6 = IQ6
                        P6 = P5 + (Q6-Q5)/T7
                        P16 = Pmin                    !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                        Q16 = Q6 + IT2 * (P16-P6) !負側の曲げ降伏点は6から16へ移動する(6から16の間が滑る部分)
                    end if
                    Gyn = P16
                case(6) !Stage 2 (荷重が負側のサイド)
                    !----負側の骨格曲線は更新なし----!
                    P4 = LD
                    Q4 = LS
                    T1 = IT1
                    P5 = P5
                    Q5 = Q5
                    T6 = T6
                    P6 = P6
                    Q6 = Q6
                    T7 = T7
                    P16 = P16
                    Q16 = Q16
                    Gyn = Gyn
                    !----正側の骨格曲線は更新あり----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    if     (LD <  IP2 .and. LD > IP5)                   then !アンカーボルトとBPLが接触していない範囲, Rule 1 or Rule 3
                        Q2 = IQ2
                        P2 = IP2
                        T2 = (Q2-Q1) / (P2-P1)
                        if     (Pmax <  IP3) then !正側の曲げ降伏非経験, Rule 1
                            P3 = IP3
                            Q3 = IQ3
                            T3 = IT3
                            P13 = P3
                            Q13 = Q3
                        else if(Pmax >= IP3) then !正側の曲げ降伏　経験, Rule 3
                            T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                            Q3 = IQ3
                            P3 = P2 + (Q3-Q2)/T3
                            P13 = Pmax                    !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                            Q13 = Q3 + IT2 * (P13-P3) !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                        end if
                    else if(LD >= IP2 .and. Pmax <  IP3) then !アンカーボルトとBPLが接触している範囲, 曲げ降伏非経験, Rule 2
                        P2 = Pmax
                        Q2 = Qmax
                        T2 = (Q2 - Q1) / (P2 - P1)
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(LD >= IP2 .and. Pmax >= IP3) then !アンカーボルトとBPLが接触している範囲, 曲げ降伏　経験, Rule 4
                        Q2 = Q1 !正側のステージ2は一時的に消失する。
                        P2 = P1 !正側のステージ2は一時的に消失する。
                        T2 = T1 !正側のステージ2は一時的に消失する。
                        P3 = Pmax
                        Q3 = IQ3
                        T3 = (Q3-Q2) / (P3-P2)
                        P13 = P3
                        Q13 = Q3
                    end if
                    Gyp = P13
                case(7) !Stage 3 (荷重が負側のサイド)
                    !----負側の骨格曲線は更新なし----!
                    P4 = LD
                    Q4 = LS
                    T1 = IT1
                    P5 = P4 !負側のStage2は一時的に消失する
                    Q5 = Q4 !負側のStage2は一時的に消失する
                    T6 = T1 !負側のStage2は一時的に消失する
                    P6 = P6
                    Q6 = Q6
                    T7 = T7
                    P16 = P16
                    Q16 = Q16
                    Gyn = Gyn
                    !----正側の骨格曲線は更新あり----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    Q2 = IQ2
                    P2 = IP2
                    T2 = (Q2-Q1)/(P2-P1)
                    if     (Pmax <  IP3) then !正側の曲げ降伏非経験, Rule 1
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(Pmax >= IP3) then !正側の曲げ降伏　経験, Rule 3
                        T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                        Q3 = IQ3
                        P3 = P2 + (Q3-Q2)/T3
                        P13 = Pmax                    !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                        Q13 = Q3 + IT2 * (P13-P3) !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                    end if
                    Gyp = P13
                case(8) !Stage 4 (荷重が負側のサイド)
                    !----負側の骨格曲線は更新なし----!
                    P4 = LD
                    Q4 = LS
                    T1 = IT1
                    P5 = P4 !負側のStage2は一時的に消失する
                    Q5 = Q4 !負側のStage2は一時的に消失する
                    T6 = T1 !負側のStage2は一時的に消失する
                    P6 = P4 !負側のStage3は一時的に消失する
                    Q6 = Q4 !負側のStage3は一時的に消失する
                    T7 = T1 !負側のStage3は一時的に消失する
                    P16 = P6
                    Q16 = Q6
                    Gyn = IP6 !負側のStage4の変形開始点は初期曲げ降伏点と同じ
                    !----正側の骨格曲線は更新あり----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    Q2 = IQ2
                    P2 = IP2
                    T2 = (Q2-Q1)/(P2-P1)
                    if     (Pmax <  IP3) then !正側の曲げ降伏非経験, Rule 1
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(Pmax >= IP3) then !正側の曲げ降伏　経験, Rule 3
                        T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                        Q3 = IQ3
                        P3 = P2 + (Q3-Q2)/T3
                        P13 = Pmax                    !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                        Q13 = Q3 + IT2 * (P13-P3) !正側の曲げ降伏点は3から13へ移動する(3から13の間が滑る部分)
                    end if
                    Gyp = P13
            end select
        end if
    end if
!--------------------------------------------------------
!Compute D from the increment
!増分変形から現在の変形を計算する。
!--------------------------------------------------------
    D = LD + DD
!--------------------------------------------------------
!Compute the current stage in the backbone curve
!骨格曲線からステージを判別する。    
!--------------------------------------------------------
    stage = 9999                        ! 9999 is the error indicator.
    if     ( D >= P13            ) then  ! 正側のStage4 (アンカーボルトの引張に伴う幾何非線形性による耐力上昇域)
        stage = 4
    else if( D >= P2            ) then  ! 正側のStage3 (アンカーボルトがBPLと接触して曲げ降伏するまでの耐力上昇域)
        stage = 3
    else if( D >= P1            ) then  ! 正側のStage2（BPLが滑り出してアンカーボルトと接触するまで)
        stage = 2
    else if( D >  P4.and.D < P1 ) then  ! Stage1 (滑り出すまで)
        stage = 1
    else if( D <= P16            ) then  ! 負側のStage4 (アンカーボルトの引張に伴う幾何非線形性による耐力上昇域)
        stage = 8
    else if( D <= P5            ) then  ! 負側のStage3 (アンカーボルトがBPLと接触して曲げ降伏するまでの耐力上昇域)
        stage = 7
    else if( D <= P4            ) then  ! 負側のStage2（BPLが滑り出してアンカーボルトと接触するまで)
        stage = 6
    else
        print *,"Stage error 9999"
        print *,"D :" ,D
        print *,"S :" ,S
        print *,"T :" ,T
        print *,"P1:",P1
        print *,"P2:",P2
        print *,"P3:",P3
        print *,"P13:",P13
        print *,"P4:",P4
        print *,"P5:",P5
        print *,"P6:",P6
        print *,"P16:",P16
        print *,"Q1:",Q1
        print *,"Q2:",Q2
        print *,"Q3:",Q3
        print *,"Q13:",Q13
        print *,"Q4:",Q4
        print *,"Q5:",Q5
        print *,"Q6:",Q6
        print *,"Q16:",Q16
        print *,"T1:",T1
        print *,"T2:",T2
        print *,"T3:",T3
        print *,"T6:",T6
        print *,"T7:",T7
        print *,"LD:" ,LD
        print *,"DD:" ,DD
        print *,"istage:",istage
        pause
        stop
    end if
!--------------------------------------------------------
!Compute next S,T based on the current backbone curve 
!骨格曲線に基いて無次元化された荷重,接線剛性比を計算する。
!--------------------------------------------------------
    select case(stage)
        case(4) ! 正側のStage4 (アンカーボルトの引張に伴う幾何非線形性による耐力上昇域)
            S = AnchorGeoNonlinear(nn,mu,IP2,Gyp,D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = IT2
        case(3) ! 正側のStage3 (アンカーボルトがBPLと接触して曲げ降伏するまでの耐力上昇域)
            if( D <  P13.and.D >= P3) then  ! 正側のStage4+Rule3 (アンカーボルトが再度引張られるまでのスリップ域)
                call LinearSlope(S,T,D,IT2,P3,Q3)
            else
                call LinearSlope(S,T,D,T3,P3,Q3)
            end if
        case(2) ! 正側のStage2（BPLが滑り出してアンカーボルトと接触するまで)
            call LinearSlope(S,T,D,T2,P2,Q2)
        case(1) ! Stage1 (滑り出すまで)
            call LinearSlope(S,T,D,T1,P1,Q1)
        case(6) ! 負側のStage2（BPLが滑り出してアンカーボルトと接触するまで)
            call LinearSlope(S,T,D,T6,P4,Q4)
        case(7) ! 負側のStage3 (アンカーボルトがBPLと接触して曲げ降伏するまでの耐力上昇域)
            if     ( D <= P6.and.D > P16 ) then  ! 負側のStage4+Rule3 (アンカーボルトが再度引張られるまでのスリップ域)
                call LinearSlope(S,T,D,IT2,P6,Q6)
            else
                call LinearSlope(S,T,D,T7,P5,Q5)
            end if
        case(8) ! 負側のStage4 (アンカーボルトの引張に伴う幾何非線形性による耐力上昇域)
            S = -1.0d0 * AnchorGeoNonlinear(nn,mu,IP2,-1.0d0*Gyn,-1.0d0*D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = 1.0d-1
    end select
!--------------------------------------------------------
!Substitute LD,LS,istage,Pmin,Pmax,Qmin,Qmax
!LD,LS,istage,Pmin,Pmax,Qmin,Qmaxを更新する。
!--------------------------------------------------------
    LD     = D
    LS     = S
    istage = stage
    if(D < Pmin) Pmin = D
    if(D > Pmax) Pmax = D
    if(S < Qmin) Qmin = S
    if(S > Qmax) Qmax = S
contains
subroutine LinearSlope(SS,TT,DDD,TTT,PPP,QQQ) !線形の骨格曲線の復元力特性計算
    implicit none
    double precision, intent( in)::DDD,TTT,PPP,QQQ
    double precision, intent(out)::SS,TT
    SS = QQQ + TTT * (DDD-PPP)
    TT = TTT
end subroutine
double precision function AnchorGeoNonlinear(nnn,mmu,Dcl,Dy,DDD,hlm,CC1,CC2,CC3,CC5) !非線形の骨格曲線の復元力特性計算(Stage4)
    implicit none
    double precision :: nnn,mmu,Dcl,Dy,DDD,hlm,CC1,CC2,CC3,CC5
    double precision :: Qs,Ns,Ms,gum,gumlrg
    gum    = datan( abs(DDD-Dcl) / abs(hlm) ) !γ
    gumlrg = datan( abs(DDD-Dy)  / abs(hlm) ) !γlrg
    Ns     = CC1 * dsin(gumlrg) * 0.25d0 * (dcos(gumlrg)**(-2.0d0))    !N/Qf={Ny^2*(h+lm)/Mp} * {sin(γlrg) / (4cos^2(γlrg))}/Qf Qfの無次元化はCCに含まれている。
    Ms     = CC2 - CC1 * (dsin(gumlrg)**2.0d0)*(dcos(gumlrg)**(-4.0d0)) * 0.0625d0 !M/(h+lm)Qf 無次元化はCCに含まれてる
    if (Ms < 0 .or. Ns >= CC5) then !塑性関係式による制約
        Ms = 0.0d0
        Ns = CC5
    end if
    Qs     = nnn*Ns*(dsin(gum)+mmu*dcos(gum)) + 2.0d0*nnn*Ms + mmu*CC3 !Q/Qf={nNsinγ+n*2M/(h+lm)+μ(nNcosγ+Pv)}/Qf 無次元化はCCに含まれてる
    AnchorGeoNonlinear = Qs
end function
double precision function FuncAlpha(Dcl,Dy,Dmax) !Rule3のStage3の剛性低減率
    implicit none
    double precision :: Dcl,Dy,Dmax
    !Dcl    :   滑り出し点
    !Dy     :   曲げ降伏点
    !Dmax   :   経験最大点
    FuncAlpha = dsqrt( abs(Dy-Dcl) / abs(Dmax-Dcl) )
    if(FuncAlpha > 1.0d0) FuncAlpha = 1.0d0
    return
end function
end subroutine