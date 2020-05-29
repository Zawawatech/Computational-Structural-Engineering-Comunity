!********************************************************
!   Subroutine  :   Driver_PlasticityModel_Uniaxial_YamashitaPinBearing
!   Purpose     :   To compute normalized material nonlinearity 
!                  (i.e. Inner force, tangent modulas)
!                   �������ʂƂ��ĕ����͓������v�Z����B
!   Reference   :   Watanabe S. and Yamashita T.,
!                   Earthquake response analysis of steel roof gymnasiums considering nonlinear
!                   restoring force characteristics of lower structures and roof bearings,
!                   Trans. AIJ ("Kibyo-shi ���\��"), AIJ, Vol.85, No.768,pp.209-218, 2020.2.
!   Remarks     :   This subroutine is based on an updated backbone curve tecnique by Yuki TERAZAWA (2020)
!                   ���̃��[�`���͍��i�Ȑ��𒀎��X�V����e�N�j�b�N��p���ď������낵�܂����B
!               :   Rule4�̃s���`���O�O�̃X���b�v��\�����邽�߂ɓ��ʂȐ܂�_13,16�������Ă��܂��B
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
!���͕ϐ���錾����
!--------------------------------------------------------
![in]
!   [integer]
!       istep   : Step number (���̃X�e�b�v�ԍ�)    
!       iter    : Iteration number in this step (���̎��ʉ񐔔ԍ�)
!   [double precision]
!       DD      : The normalized incremental deformation of n-th iteration in this step (����̖������������ό`)
!       IP1     : (Positive deformation) Break point 1 in the initial backbone curve (�����̍��i�Ȑ��܂�_1,�����̊���o��,�ό`)
!       IP2     : (Positive deformation) Break point 2 in the initial backbone curve (�����̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g�ڐG,�ό`)
!       IP3     : (Positive deformation) Break point 3 in the initial backbone curve (�����̍��i�Ȑ��܂�_3,�����̃A���J�[�{���g�Ȃ��~��,�ό`)
!       IP4     : (Negative deformation) Break point 4 in the initial backbone curve (�����̍��i�Ȑ��܂�_4,�����̊���o��,�ό`)
!       IP5     : (Negative deformation) Break point 5 in the initial backbone curve (�����̍��i�Ȑ��܂�_5,������BPL�ƃA���J�[�{���g�ڐG,�ό`)
!       IP6     : (Negative deformation) Break point 6 in the initial backbone curve (�����̍��i�Ȑ��܂�_6,�����̃A���J�[�{���g�Ȃ��~��,�ό`)    
!       IQ1     : (Positive force) Break point 1 in the initial backbone curve (�����̍��i�Ȑ��܂�_1,�����̊���o��,������)
!       IQ2     : (Positive force) Break point 2 in the initial backbone curve (�����̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g�ڐG,������)
!       IQ3     : (Positive force) Break point 3 in the initial backbone curve (�����̍��i�Ȑ��܂�_3,�����̃A���J�[�{���g�Ȃ��~��,������)
!       IQ4     : (Negative force) Break point 4 in the initial backbone curve (�����̍��i�Ȑ��܂�_4,�����̊���o��,������)
!       IQ5     : (Negative force) Break point 5 in the initial backbone curve (�����̍��i�Ȑ��܂�_5,������BPL�ƃA���J�[�{���g�ڐG,������)
!       IQ6     : (Negative force) Break point 6 in the initial backbone curve (�����̍��i�Ȑ��܂�_6,�����̃A���J�[�{���g�Ȃ��~��,������)        
!       IT1     : (Tangent modulas) Slope 1 in the initial backbone curve (�����̍��i�Ȑ��̌��z1,�R�����f����Stage 1 )
!       IT2     : (Tangent modulas) Slope 2 in the initial backbone curve (�����̍��i�Ȑ��̌��z2,�R�����f����Stage 2 )
!       IT3     : (Tangent modulas) Slope 3 in the initial backbone curve (�����̍��i�Ȑ��̌��z3,�R�����f����Stage 3 )
!       nn      : the number of anchorbolt. (�A���J�[�{���g�{��)
!       mu      : friction coeficient (���C�W��)
!       C1      : C1 = {Ny*Ny*(h+lm)/Mp}/Qf = N/Qf
!       C2      : C2 = {Mp/(h+lm)}/Qf
!       C3      : C3 = Pv/Qf
!       C4      : C4 = (h+lm)/��f
!       C5      : C5 = Ny/Qf
![inout]
!   [integer]    
!       istage  : ((n-1)-th iteration) The last stage  (�O��̃X�e�[�W�ԍ�)    
!   [double precision]
!       LD      : ((n-1)-th iteration) The last normalized deformation    (�O��̖��������ό`)
!       LS      : ((n-1)-th iteration) The last normalized stress         (�O��̖��������׏d)
!       D       : (n-th iteration)     The current normalized deformation (����̖��������ό`)
!       S       : (n-th iteration)     The current normalized stress      (����̖��������׏d)
!       T       : (n-th iteration)     The current tangent modulas        (����̐ڐ�������)
!       P1      : (Positive deformation) Break point 1 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_1,�����̊���o��,�ό`)
!       P2      : (Positive deformation) Break point 2 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g�ڐG,�ό`)
!       P3      : (Positive deformation) Break point 3 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_3,�����̃A���J�[�{���g�Ȃ��~��,�ό`)
!       P13     : (Positive deformation) Break point13 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_13,�čډ׎��̐����̃A���J�[�{���g�̃s���`���O�J�n�_,�ό`)
!       P4      : (Negative deformation) Break point 4 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_4,�����̊���o��,�ό`)
!       P5      : (Negative deformation) Break point 5 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_5,������BPL�ƃA���J�[�{���g�ڐG,�ό`)
!       P6      : (Negative deformation) Break point 6 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_6,�����̃A���J�[�{���g�Ȃ��~��,�ό`)
!       P16     : (Negative deformation) Break point16 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_16,�čډ׎��̐����̃A���J�[�{���g�̃s���`���O�J�n�_,�ό`)
!       Q1      : (Positive force) Break point 1 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_1,�����̊���o��,������)
!       Q2      : (Positive force) Break point 2 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g�ڐG,������)
!       Q3      : (Positive force) Break point 3 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_3,�����̃A���J�[�{���g�Ȃ��~��,������)
!       Q13     : (Positive force) Break point13 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_13,�čډ׎��̐����̃A���J�[�{���g�̃s���`���O�J�n�_,������)
!       Q4      : (Negative force) Break point 4 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_4,�����̊���o��,������)
!       Q5      : (Negative force) Break point 5 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_5,������BPL�ƃA���J�[�{���g�ڐG,������)
!       Q6      : (Negative force) Break point 6 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_6,�����̃A���J�[�{���g�Ȃ��~��,������)
!       Q16     : (Negative force) Break point16 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_16,�čډ׎��̐����̃A���J�[�{���g�̃s���`���O�J�n�_,������)   
!       T1      : (Tangent modulas) Slope 1 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z1,�R�����f����Stage 1,����,����O)
!       T2      : (Tangent modulas) Slope 2 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z2,�R�����f����Stage 2,����,���蒆)
!       T3      : (Tangent modulas) Slope 3 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z3,�R�����f����Stage 3,����,�A���J�[�Ȃ��ό`��)
!       T6      : (Tangent modulas) Slope 6 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z6,�R�����f����Stage 2,����,���蒆)
!       T7      : (Tangent modulas) Slope 7 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z7,�R�����f����Stage 3,����,�A���J�[�Ȃ��ό`��)
!       Pmax    : Maximum value of experienced deformation (�o���ő�ό`)
!       Pmin    : Minimum value of experienced deformation (�o���ŏ��ό`)    
!       Qmax    : Maximum value of experienced force (�o���ő啜����)
!       Qmin    : Minimum value of experienced force (�o���ŏ�������)
!       Gyp     : �s���`���O���n�܂�ό`�_(����)
!       Gyn     : �s���`���O���n�܂�ό`�_(����)
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
!���[�J���ϐ��̐錾    
!--------------------------------------------------------
![integer]
!   stage   : Stage number in the backbone curve (���i�Ȑ����̃X�e�[�W�敪)
!--------------------------------------------------------    
    integer             :: stage
!--------------------------------------------------------
!Update backbone curve if iter = 0 
!����0��ڂȂ獜�i�Ȑ����X�V����B
!--------------------------------------------------------
    if(iter == 0) then
        if (istep == 1) then    !��͊J�n���_�̏����l���
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
            select case(istage) !�O��̃X�e�b�v�̍ŏI���ʎ��_�̉׏d�ό`�֌W�̍��W���獡��̃X�e�b�v�̍��i�Ȑ������B
                case(2) !Stage 2 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
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
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    if     (LD <  IP2 .and. LD > IP5)                   then !�A���J�[�{���g��BPL���ڐG���Ă��Ȃ��͈�, Rule 1 or Rule 3
                        Q5 = IQ5
                        P5 = IP5
                        T6 = (Q4-Q5) / (P4-P5)
                        if     (Pmin >  IP6) then !�����̋Ȃ��~����o��, Rule 1
                            P6 = IP6
                            Q6 = IQ6
                            T7 = IT3
                            P16 = P6
                            Q16 = Q6
                        else if(Pmin <= IP6) then !�����̋Ȃ��~���@�o��, Rule 3
                            T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                            Q6 = IQ6
                            P6 = P5 + (Q6-Q5)/T7
                            P16 = Pmin                    !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                            Q16 = Q6 + IT2 * (P16-P6) !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                        end if
                    else if(LD <= IP5 .and. Pmin >  IP6) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~����o��, Rule 2
                        P5 = Pmin
                        Q5 = Qmin
                        T6 = (Q4 - Q5) / (P4 - P5)
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(LD <= IP5 .and. Pmin <= IP6) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~���@�o��, Rule 4
                        Q5 = Q4 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        P5 = P4 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        T6 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        P6 = Pmin
                        Q6 = IQ6
                        T7 = (Q5-Q6) / (P5-P6)
                        P16 = P6
                        Q16 = Q6
                    end if
                    Gyn = P16
                case(3) !Stage 3 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !������Stage2�͈ꎞ�I�ɏ�������
                    Q2 = Q1 !������Stage2�͈ꎞ�I�ɏ�������
                    T2 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P3 = P3
                    Q3 = Q3
                    T3 = T3
                    P13 = P13
                    Q13 = Q13
                    Gyp = Gyp
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    Q5 = IQ5
                    P5 = IP5
                    T6 = (Q4-Q5)/(P4-P5)
                    if     (Pmin >  IP6) then !�����̋Ȃ��~����o��, Rule 1
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(Pmin <= IP6) then !�����̋Ȃ��~���@�o��, Rule 3
                        T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                        Q6 = IQ6
                        P6 = P5 + (Q6-Q5)/T7
                        P16 = Pmin                    !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                        Q16 = Q6 + IT2 * (P16-P6) !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                    end if
                    Gyn = P16
                case(4) !Stage 4 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !������Stage2�͈ꎞ�I�ɏ�������
                    Q2 = Q1 !������Stage2�͈ꎞ�I�ɏ�������
                    T2 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P3 = P1 !������Stage3�͈ꎞ�I�ɏ�������
                    Q3 = Q1 !������Stage3�͈ꎞ�I�ɏ�������
                    T3 = T1 !������Stage3�͈ꎞ�I�ɏ�������
                    P13 = P3
                    Q13 = Q3
                    Gyp = IP3 !������Stage4�̕ό`�J�n�_�͏����Ȃ��~���_�Ɠ���
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q4 = IQ4
                    P4 = P1 + (Q4-Q1)/T1
                    Q5 = IQ5
                    P5 = IP5
                    T6 = (Q4-Q5)/(P4-P5)
                    if     (Pmin >  IP6) then !�����̋Ȃ��~����o��, Rule 1
                        P6 = IP6
                        Q6 = IQ6
                        T7 = IT3
                        P16 = P6
                        Q16 = Q6
                    else if(Pmin <= IP6) then !�����̋Ȃ��~���@�o��, Rule 3
                        T7 = FuncAlpha(IP2,IP3,-1.0d0*Pmin) * IT3
                        Q6 = IQ6
                        P6 = P5 + (Q6-Q5)/T7
                        P16 = Pmin                    !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                        Q16 = Q6 + IT2 * (P16-P6) !�����̋Ȃ��~���_��6����16�ֈړ�����(6����16�̊Ԃ����镔��)
                    end if
                    Gyn = P16
                case(6) !Stage 2 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
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
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    if     (LD <  IP2 .and. LD > IP5)                   then !�A���J�[�{���g��BPL���ڐG���Ă��Ȃ��͈�, Rule 1 or Rule 3
                        Q2 = IQ2
                        P2 = IP2
                        T2 = (Q2-Q1) / (P2-P1)
                        if     (Pmax <  IP3) then !�����̋Ȃ��~����o��, Rule 1
                            P3 = IP3
                            Q3 = IQ3
                            T3 = IT3
                            P13 = P3
                            Q13 = Q3
                        else if(Pmax >= IP3) then !�����̋Ȃ��~���@�o��, Rule 3
                            T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                            Q3 = IQ3
                            P3 = P2 + (Q3-Q2)/T3
                            P13 = Pmax                    !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                            Q13 = Q3 + IT2 * (P13-P3) !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                        end if
                    else if(LD >= IP2 .and. Pmax <  IP3) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~����o��, Rule 2
                        P2 = Pmax
                        Q2 = Qmax
                        T2 = (Q2 - Q1) / (P2 - P1)
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(LD >= IP2 .and. Pmax >= IP3) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~���@�o��, Rule 4
                        Q2 = Q1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        P2 = P1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        T2 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        P3 = Pmax
                        Q3 = IQ3
                        T3 = (Q3-Q2) / (P3-P2)
                        P13 = P3
                        Q13 = Q3
                    end if
                    Gyp = P13
                case(7) !Stage 3 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P4 = LD
                    Q4 = LS
                    T1 = IT1
                    P5 = P4 !������Stage2�͈ꎞ�I�ɏ�������
                    Q5 = Q4 !������Stage2�͈ꎞ�I�ɏ�������
                    T6 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P6 = P6
                    Q6 = Q6
                    T7 = T7
                    P16 = P16
                    Q16 = Q16
                    Gyn = Gyn
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    Q2 = IQ2
                    P2 = IP2
                    T2 = (Q2-Q1)/(P2-P1)
                    if     (Pmax <  IP3) then !�����̋Ȃ��~����o��, Rule 1
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(Pmax >= IP3) then !�����̋Ȃ��~���@�o��, Rule 3
                        T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                        Q3 = IQ3
                        P3 = P2 + (Q3-Q2)/T3
                        P13 = Pmax                    !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                        Q13 = Q3 + IT2 * (P13-P3) !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                    end if
                    Gyp = P13
                case(8) !Stage 4 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P4 = LD
                    Q4 = LS
                    T1 = IT1
                    P5 = P4 !������Stage2�͈ꎞ�I�ɏ�������
                    Q5 = Q4 !������Stage2�͈ꎞ�I�ɏ�������
                    T6 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P6 = P4 !������Stage3�͈ꎞ�I�ɏ�������
                    Q6 = Q4 !������Stage3�͈ꎞ�I�ɏ�������
                    T7 = T1 !������Stage3�͈ꎞ�I�ɏ�������
                    P16 = P6
                    Q16 = Q6
                    Gyn = IP6 !������Stage4�̕ό`�J�n�_�͏����Ȃ��~���_�Ɠ���
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P4 + (Q1-Q4)/T1
                    Q2 = IQ2
                    P2 = IP2
                    T2 = (Q2-Q1)/(P2-P1)
                    if     (Pmax <  IP3) then !�����̋Ȃ��~����o��, Rule 1
                        P3 = IP3
                        Q3 = IQ3
                        T3 = IT3
                        P13 = P3
                        Q13 = Q3
                    else if(Pmax >= IP3) then !�����̋Ȃ��~���@�o��, Rule 3
                        T3 = FuncAlpha(IP2,IP3,Pmax) * IT3
                        Q3 = IQ3
                        P3 = P2 + (Q3-Q2)/T3
                        P13 = Pmax                    !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                        Q13 = Q3 + IT2 * (P13-P3) !�����̋Ȃ��~���_��3����13�ֈړ�����(3����13�̊Ԃ����镔��)
                    end if
                    Gyp = P13
            end select
        end if
    end if
!--------------------------------------------------------
!Compute D from the increment
!�����ό`���猻�݂̕ό`���v�Z����B
!--------------------------------------------------------
    D = LD + DD
!--------------------------------------------------------
!Compute the current stage in the backbone curve
!���i�Ȑ�����X�e�[�W�𔻕ʂ���B    
!--------------------------------------------------------
    stage = 9999                        ! 9999 is the error indicator.
    if     ( D >= P13            ) then  ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
        stage = 4
    else if( D >= P2            ) then  ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
        stage = 3
    else if( D >= P1            ) then  ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
        stage = 2
    else if( D >  P4.and.D < P1 ) then  ! Stage1 (����o���܂�)
        stage = 1
    else if( D <= P16            ) then  ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
        stage = 8
    else if( D <= P5            ) then  ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
        stage = 7
    else if( D <= P4            ) then  ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
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
!���i�Ȑ��Ɋ�Ė����������ꂽ�׏d,�ڐ���������v�Z����B
!--------------------------------------------------------
    select case(stage)
        case(4) ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
            S = AnchorGeoNonlinear(nn,mu,IP2,Gyp,D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = IT2
        case(3) ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
            if( D <  P13.and.D >= P3) then  ! ������Stage4+Rule3 (�A���J�[�{���g���ēx��������܂ł̃X���b�v��)
                call LinearSlope(S,T,D,IT2,P3,Q3)
            else
                call LinearSlope(S,T,D,T3,P3,Q3)
            end if
        case(2) ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
            call LinearSlope(S,T,D,T2,P2,Q2)
        case(1) ! Stage1 (����o���܂�)
            call LinearSlope(S,T,D,T1,P1,Q1)
        case(6) ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
            call LinearSlope(S,T,D,T6,P4,Q4)
        case(7) ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
            if     ( D <= P6.and.D > P16 ) then  ! ������Stage4+Rule3 (�A���J�[�{���g���ēx��������܂ł̃X���b�v��)
                call LinearSlope(S,T,D,IT2,P6,Q6)
            else
                call LinearSlope(S,T,D,T7,P5,Q5)
            end if
        case(8) ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
            S = -1.0d0 * AnchorGeoNonlinear(nn,mu,IP2,-1.0d0*Gyn,-1.0d0*D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = 1.0d-1
    end select
!--------------------------------------------------------
!Substitute LD,LS,istage,Pmin,Pmax,Qmin,Qmax
!LD,LS,istage,Pmin,Pmax,Qmin,Qmax���X�V����B
!--------------------------------------------------------
    LD     = D
    LS     = S
    istage = stage
    if(D < Pmin) Pmin = D
    if(D > Pmax) Pmax = D
    if(S < Qmin) Qmin = S
    if(S > Qmax) Qmax = S
contains
subroutine LinearSlope(SS,TT,DDD,TTT,PPP,QQQ) !���`�̍��i�Ȑ��̕����͓����v�Z
    implicit none
    double precision, intent( in)::DDD,TTT,PPP,QQQ
    double precision, intent(out)::SS,TT
    SS = QQQ + TTT * (DDD-PPP)
    TT = TTT
end subroutine
double precision function AnchorGeoNonlinear(nnn,mmu,Dcl,Dy,DDD,hlm,CC1,CC2,CC3,CC5) !����`�̍��i�Ȑ��̕����͓����v�Z(Stage4)
    implicit none
    double precision :: nnn,mmu,Dcl,Dy,DDD,hlm,CC1,CC2,CC3,CC5
    double precision :: Qs,Ns,Ms,gum,gumlrg
    gum    = datan( abs(DDD-Dcl) / abs(hlm) ) !��
    gumlrg = datan( abs(DDD-Dy)  / abs(hlm) ) !��lrg
    Ns     = CC1 * dsin(gumlrg) * 0.25d0 * (dcos(gumlrg)**(-2.0d0))    !N/Qf={Ny^2*(h+lm)/Mp} * {sin(��lrg) / (4cos^2(��lrg))}/Qf Qf�̖���������CC�Ɋ܂܂�Ă���B
    Ms     = CC2 - CC1 * (dsin(gumlrg)**2.0d0)*(dcos(gumlrg)**(-4.0d0)) * 0.0625d0 !M/(h+lm)Qf ����������CC�Ɋ܂܂�Ă�
    if (Ms < 0 .or. Ns >= CC5) then !�Y���֌W���ɂ�鐧��
        Ms = 0.0d0
        Ns = CC5
    end if
    Qs     = nnn*Ns*(dsin(gum)+mmu*dcos(gum)) + 2.0d0*nnn*Ms + mmu*CC3 !Q/Qf={nNsin��+n*2M/(h+lm)+��(nNcos��+Pv)}/Qf ����������CC�Ɋ܂܂�Ă�
    AnchorGeoNonlinear = Qs
end function
double precision function FuncAlpha(Dcl,Dy,Dmax) !Rule3��Stage3�̍����ጸ��
    implicit none
    double precision :: Dcl,Dy,Dmax
    !Dcl    :   ����o���_
    !Dy     :   �Ȃ��~���_
    !Dmax   :   �o���ő�_
    FuncAlpha = dsqrt( abs(Dy-Dcl) / abs(Dmax-Dcl) )
    if(FuncAlpha > 1.0d0) FuncAlpha = 1.0d0
    return
end function
end subroutine