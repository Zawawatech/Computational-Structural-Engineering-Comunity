!********************************************************
!   Subroutine  :   Driver_PlasticityModel_Uniaxial_YamashitaPinRollerBearing
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
subroutine Driver_PlasticityModel_Uniaxial_YamashitaPinRollerBearing(&
            istep,iter,istage,&
            LD,LS,&
            D,S,T,&
            DD,&
            P1,P2,P3,P4,P14,P5,P6,P7,P8,P18,&
            Q1,Q2,Q3,Q4,Q14,Q5,Q6,Q7,Q8,Q18,&
            T1,T2,T3,T4,T6,T7,T8,&
            IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8,&
            IQ1,IQ2,IQ3,IQ4,IQ5,IQ6,IQ7,IQ8,&
            IT1,IT2,IT3,IT4,&
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
!       IP1     : (Positive deformation) Break point 1 in the initial backbone curve (�����̍��i�Ȑ��܂�_1,������BPL����o��,�ό`)
!       IP2     : (Positive deformation) Break point 2 in the initial backbone curve (�����̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,�ό`)    
!       IP3     : (Positive deformation) Break point 3 in the initial backbone curve (�����̍��i�Ȑ��܂�_3,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,�ό`)
!       IP4     : (Positive deformation) Break point 4 in the initial backbone curve (�����̍��i�Ȑ��܂�_4,�����̃A���J�[�{���g�Ȃ��~��,�ό`)
!       IP5     : (Negative deformation) Break point 5 in the initial backbone curve (�����̍��i�Ȑ��܂�_5,������BPL����o��,�ό`)
!       IP6     : (Negative deformation) Break point 6 in the initial backbone curve (�����̍��i�Ȑ��܂�_6,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,�ό`)        
!       IP7     : (Negative deformation) Break point 7 in the initial backbone curve (�����̍��i�Ȑ��܂�_7,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,�ό`)
!       IP8     : (Negative deformation) Break point 8 in the initial backbone curve (�����̍��i�Ȑ��܂�_8,�����̃A���J�[�{���g�Ȃ��~��,�ό`)    
!       IQ1     : (Positive force) Break point 1 in the initial backbone curve (�����̍��i�Ȑ��܂�_1,������BPL����o��,������)
!       IQ2     : (Positive force) Break point 2 in the initial backbone curve (�����̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,������)        
!       IQ3     : (Positive force) Break point 3 in the initial backbone curve (�����̍��i�Ȑ��܂�_3,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,������)
!       IQ4     : (Positive force) Break point 4 in the initial backbone curve (�����̍��i�Ȑ��܂�_4,�����̃A���J�[�{���g�Ȃ��~��,������)
!       IQ5     : (Negative force) Break point 5 in the initial backbone curve (�����̍��i�Ȑ��܂�_5,������BPL����o��,������)
!       IQ6     : (Negative force) Break point 6 in the initial backbone curve (�����̍��i�Ȑ��܂�_6,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,������)
!       IQ7     : (Negative force) Break point 7 in the initial backbone curve (�����̍��i�Ȑ��܂�_7,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,������)
!       IQ8     : (Negative force) Break point 8 in the initial backbone curve (�����̍��i�Ȑ��܂�_8,�����̃A���J�[�{���g�Ȃ��~��,������)        
!       IT1     : (Tangent modulas) Slope 1 in the initial backbone curve (�����̍��i�Ȑ��̌��z1,�R�����f����Stage 1 )
!       IT2     : (Tangent modulas) Slope 2 in the initial backbone curve (�����̍��i�Ȑ��̌��z2,�R�����f����Stage 2 )
!       IT3     : (Tangent modulas) Slope 3 in the initial backbone curve (�����̍��i�Ȑ��̌��z3,�R�����f����Stage 2')
!       IT4     : (Tangent modulas) Slope 4 in the initial backbone curve (�����̍��i�Ȑ��̌��z4,�R�����f����Stage 3 )
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
!       P1      : (Positive deformation) Break point 1 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_1,������BPL����o��,�ό`)
!       P2      : (Positive deformation) Break point 2 in the curent backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,�ό`)    
!       P3      : (Positive deformation) Break point 3 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_3,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,�ό`)
!       P4      : (Positive deformation) Break point 4 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_4,�����̃A���J�[�{���g�Ȃ��~��,�ό`)
!       P14     : (Positive deformation) Break point14 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_14,�����̃A���J�[�{���g�Ȃ��~��,�ό`)    
!       P5      : (Negative deformation) Break point 5 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_5,������BPL����o��,�ό`)
!       P6      : (Negative deformation) Break point 6 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_6,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,�ό`)        
!       P7      : (Negative deformation) Break point 7 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_7,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,�ό`)
!       P8      : (Negative deformation) Break point 8 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_8,�����̃A���J�[�{���g�Ȃ��~��,�ό`)    
!       P18     : (Negative deformation) Break point18 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_18,�����̃A���J�[�{���g�Ȃ��~��,�ό`)        
!       Q1      : (Positive force) Break point 1 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_1,������BPL����o��,������)
!       Q2      : (Positive force) Break point 2 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_2,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,������)        
!       Q3      : (Positive force) Break point 3 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_3,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,������)
!       Q4      : (Positive force) Break point 4 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_4,�����̃A���J�[�{���g�Ȃ��~��,������)
!       Q14     : (Positive force) Break point14 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_14,�����̃A���J�[�{���g�Ȃ��~��,������)    
!       Q5      : (Negative force) Break point 5 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_5,������BPL����o��,������)
!       Q6      : (Negative force) Break point 6 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_6,������BPL�ƃA���J�[�{���g���ڐG��ɕ~�v���[�g�ƕ~�����^�������C�𐶂���,������)
!       Q7      : (Negative force) Break point 7 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_7,�����̕~�v���[�g�ƕ~�����^��������o����BPL�ƃA���J�[�{���g�ڐG�ɂ��Ȃ��͂𔭊�,������)
!       Q8      : (Negative force) Break point 8 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_8,�����̃A���J�[�{���g�Ȃ��~��,������)        
!       Q18     : (Negative force) Break point18 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��܂�_18,�����̃A���J�[�{���g�Ȃ��~��,������)            
!       T1      : (Tangent modulas) Slope 1 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z1,�R�����f����Stage 1,����,����O)
!       T2      : (Tangent modulas) Slope 2 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z2,�R�����f����Stage 2,����,���蒆)
!       T3      : (Tangent modulas) Slope 3 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z3,�R�����f����Stage 2',����,�~�v���[�g�ƃ����^�����C��)    
!       T4      : (Tangent modulas) Slope 4 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z4,�R�����f����Stage 3,����,�A���J�[�Ȃ��ό`��)
!       T6      : (Tangent modulas) Slope 6 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z6,�R�����f����Stage 2,����,���蒆)
!       T7      : (Tangent modulas) Slope 7 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z7,�R�����f����Stage 2',����,�~�v���[�g�ƃ����^�����C��)        
!       T8      : (Tangent modulas) Slope 8 in the current backbone curve (���̃X�e�b�v�̍��i�Ȑ��̌��z8,�R�����f����Stage 3,����,�A���J�[�Ȃ��ό`��)
!       Pmax    : Maximum value of experienced deformation (�o���ő�ό`)
!       Pmin    : Minimum value of experienced deformation (�o���ŏ��ό`)    
!       Qmax    : Maximum value of experienced force (�o���ő啜����)
!       Qmin    : Minimum value of experienced force (�o���ŏ�������)
!       Gyp     : �s���`���O���n�܂�ό`�_(����)
!       Gyn     : �s���`���O���n�܂�ό`�_(����)
!--------------------------------------------------------
    integer,            intent(   in)   :: iter,istep
    double precision,   intent(   in)   :: DD
    double precision,   intent(   in)   :: IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8
    double precision,   intent(   in)   :: IQ1,IQ2,IQ3,IQ4,IQ5,IQ6,IQ7,IQ8
    double precision,   intent(   in)   :: IT1,IT2,IT3,IT4
    double precision,   intent(   in)   :: nn,mu,C1,C2,C3,C4,C5
    integer,            intent(inout)   :: istage
    double precision,   intent(inout)   :: LD,LS
    double precision,   intent(inout)   :: D,S,T
    double precision,   intent(inout)   :: P1,P2,P3,P4,P14,P5,P6,P7,P8,P18
    double precision,   intent(inout)   :: Q1,Q2,Q3,Q4,Q14,Q5,Q6,Q7,Q8,Q18
    double precision,   intent(inout)   :: T1,T2,T3,T4,T6,T7,T8
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
            P4 = IP4
            P14 = IP4
            P5 = IP5
            P6 = IP6
            P7 = IP7
            P8 = IP8
            P18 = IP8
            Q1 = IQ1
            Q2 = IQ2
            Q3 = IQ3
            Q4 = IQ4
            Q14 = IQ4
            Q5 = IQ5
            Q6 = IQ6
            Q7 = IQ7
            Q8 = IQ8
            Q18 = IQ8
            T1 = IT1
            T2 = IT2
            T3 = IT3
            T4 = IT4
            T6 = IT2
            T7 = IT3
            T8 = IT4
            Pmin = 0.0d0
            Pmax = 0.0d0
            Qmin = 0.0d0
            Qmax = 0.0d0
            Gyn = IP8
            Gyp = IP4
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
                    P4 = P4
                    Q4 = Q4
                    T4 = T4
                    P14 = P14
                    Q14 = Q14
                    Gyp = Gyp
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q5 = IQ5
                    P5 = P1 + (Q5-Q1)/T1
                    if     (LD <  IP3 .and. LD > IP7)                   then !�A���J�[�{���g��BPL���ڐG���Ă��Ȃ��͈�, Rule 1 or Rule 3
                        P6 = IP6
                        Q6 = IQ6
                        T6 = (Q5-Q6) / (P5-P6)
                        Q7 = IQ7
                        P7 = IP7
                        T7 = IT3
                        if     (Pmin >  IP8) then !�����̋Ȃ��~����o��, Rule 1
                            P8 = IP8
                            Q8 = IQ8
                            T8 = IT4
                            P18 = P8
                            Q18 = Q8
                        else if(Pmin <= IP8) then !�����̋Ȃ��~���@�o��, Rule 3
                            T8 = FuncAlpha(IP2,IP4,-1.0d0*Pmin) * IT4
                            Q8 = IQ8
                            P8 = P7 + (Q8-Q7)/T8
                            P18 = Pmin                   !�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                            Q18 = Q8 + IT2 * (P18-P8)!�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                        end if
                    else if(LD <= IP7 .and. Pmin >  IP8) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~����o��, Rule 2
                        P6 = Pmin
                        Q6 = Qmin
                        T6 = (Q5-Q6) / (P5-P6)
                        P7 = P6
                        Q7 = Q6
                        T7 = T6
                        P8 = IP8
                        Q8 = IQ8
                        T8 = IT4
                        P18 = P8
                        Q18 = Q8
                    else if(LD <= IP7 .and. Pmin <= IP8) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~���@�o��, Rule 4
                        Q6 = Q5 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        P6 = P5 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        T6 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������
                        Q7 = Q5 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������
                        P7 = P5 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������
                        T7 = T1 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������
                        P8 = Pmin
                        Q8 = IQ8
                        T8 = (Q7-Q8) / (P7-P8)
                        P18 = P8
                        Q18 = Q8
                    end if
                    Gyn = P18
                case(3) !Stage 2' (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    Q2 = Q1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    T2 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    P3 = P3
                    Q3 = Q3
                    T3 = T3
                    P4 = P4
                    Q4 = Q4
                    T4 = T4
                    P14 = P14
                    Q14 = Q14
                    Gyp = Gyp
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q5 = IQ5
                    P5 = P1 + (Q5-Q1)/T1
                    P6 = IP6
                    Q6 = IQ6
                    T6 = (Q5-Q6)/(P5-P6)
                    Q7 = IQ7
                    P7 = IP7
                    T7 = IT3
                    if     (Pmin >  IP8) then !�����̋Ȃ��~����o��, Rule 1
                        P8 = IP8
                        Q8 = IQ8
                        T8 = IT4
                        P18 = P8
                        Q18 = Q8
                    else if(Pmin <= IP8) then !�����̋Ȃ��~���@�o��, Rule 3
                        T8 = FuncAlpha(IP2,IP4,-1.0d0*Pmin) * IT4
                        Q8 = IQ8
                        P8 = P7 + (Q8-Q7)/T8
                        P18 = Pmin                   !�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                        Q18 = Q8 + IT2 * (P18-P8)!�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                    end if
                    Gyn = P18
                case(4) !Stage 3 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !������Stage2�͈ꎞ�I�ɏ�������
                    Q2 = Q1 !������Stage2�͈ꎞ�I�ɏ�������
                    T2 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P3 = P1 !������Stage2'�͈ꎞ�I�ɏ�������
                    Q3 = Q1 !������Stage2'�͈ꎞ�I�ɏ�������
                    T3 = T1 !������Stage2'�͈ꎞ�I�ɏ�������
                    P4 = P4
                    Q4 = Q4
                    T4 = T4
                    P14 = P14
                    Q14 = Q14
                    Gyp = Gyp
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q5 = IQ5
                    P5 = P1 + (Q5-Q1)/T1
                    P6 = IP6
                    Q6 = IQ6
                    T6 = (Q5-Q6)/(P5-P6)
                    Q7 = IQ7
                    P7 = IP7
                    T7 = IT3
                    if     (Pmin >  IP8) then !�����̋Ȃ��~����o��, Rule 1
                        P8 = IP8
                        Q8 = IQ8
                        T8 = IT4
                        P18 = P8
                        Q18 = Q8
                    else if(Pmin <= IP8) then !�����̋Ȃ��~���@�o��, Rule 3
                        T8 = FuncAlpha(IP2,IP4,-1.0d0*Pmin) * IT4
                        Q8 = IQ8
                        P8 = P7 + (Q8-Q7)/T8
                        P18 = Pmin                   !�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                        Q18 = Q8 + IT2 * (P18-P8)!�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                    end if
                    Gyn = P18
                case(5) !Stage 4 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P1 = LD
                    Q1 = LS
                    T1 = IT1
                    P2 = P1 !������Stage2�͈ꎞ�I�ɏ�������
                    Q2 = Q1 !������Stage2�͈ꎞ�I�ɏ�������
                    T2 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P3 = P1 !������Stage2'�͈ꎞ�I�ɏ�������
                    Q3 = Q1 !������Stage2'�͈ꎞ�I�ɏ�������
                    T3 = T1 !������Stage2'�͈ꎞ�I�ɏ�������
                    P4 = P1 !������Stage3�͈ꎞ�I�ɏ�������
                    Q4 = Q1 !������Stage3�͈ꎞ�I�ɏ�������
                    T4 = T1 !������Stage3�͈ꎞ�I�ɏ�������
                    P14 = P4
                    Q14 = Q4
                    Gyp = IP4 !������Stage4�̕ό`�J�n�_�͏����Ȃ��~���_�Ɠ���
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q5 = IQ5
                    P5 = P1 + (Q5-Q1)/T1
                    P6 = IP6
                    Q6 = IQ6
                    T6 = (Q5-Q6)/(P5-P6)
                    Q7 = IQ7
                    P7 = IP7
                    T7 = IT3
                    if     (Pmin >  IP8) then !�����̋Ȃ��~����o��, Rule 1
                        P8 = IP8
                        Q8 = IQ8
                        T8 = IT4
                        P18 = P8
                        Q18 = Q8
                    else if(Pmin <= IP8) then !�����̋Ȃ��~���@�o��, Rule 3
                        T8 = FuncAlpha(IP2,IP4,-1.0d0*Pmin) * IT4
                        Q8 = IQ8
                        P8 = P7 + (Q8-Q7)/T8
                        P18 = Pmin                   !�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                        Q18 = Q8 + IT2 * (P18-P8)!�����̋Ȃ��~���_��8����18�ֈړ�����(8����18�̊Ԃ����镔��)
                    end if
                    Gyn = P18
                case(6) !Stage 2 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P5 = LD
                    Q5 = LS
                    T1 = IT1
                    P6 = P6
                    Q6 = Q6
                    T6 = T6
                    P7 = P7
                    Q7 = Q7
                    T7 = T7
                    P8 = P8
                    Q8 = Q8
                    T8 = T8
                    P18 = P18
                    Q18 = Q18
                    Gyn = Gyn
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P5 + (Q1-Q5)/T1
                    if     (LD <  IP3 .and. LD > IP7)                   then !�A���J�[�{���g��BPL���ڐG���Ă��Ȃ��͈�, Rule 1 or Rule 3
                        P2 = IP2
                        Q2 = IQ2
                        T2 = (Q2-Q1) / (P2-P1)
                        Q3 = IQ3
                        P3 = IP3
                        T3 = IT3
                        if     (Pmax <  IP4) then !�����̋Ȃ��~����o��, Rule 1
                            P4 = IP4
                            Q4 = IQ4
                            T4 = IT4
                            P14 = P4
                            Q14 = Q4
                        else if(Pmax >= IP4) then !�����̋Ȃ��~���@�o��, Rule 3
                            T4 = FuncAlpha(IP2,IP4,Pmax) * IT4
                            Q4 = IQ4
                            P4 = P3 + (Q4-Q3)/T4
                            P14 = Pmax                    !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                            Q14 = Q4 + IT2 * (P14-P4) !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                        end if
                    else if(LD >= IP3 .and. Pmax <  IP4) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~����o��, Rule 2
                        P2 = Pmax
                        Q2 = Qmax
                        T2 = (Q2-Q1) / (P2-P1)
                        P3 = P2
                        Q3 = Q2
                        T3 = T2
                        P4 = IP4
                        Q4 = IQ4
                        T4 = IT4
                        P14 = P4
                        Q14 = Q4
                    else if(LD >= IP3 .and. Pmax >= IP4) then !�A���J�[�{���g��BPL���ڐG���Ă���͈�, �Ȃ��~���@�o��, Rule 4
                        P2 = P1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        Q2 = Q1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        T2 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                        Q3 = Q1 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������B
                        P3 = P1 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������B
                        T3 = T1 !�����̃X�e�[�W2'�͈ꎞ�I�ɏ�������B
                        P4 = Pmax
                        Q4 = IQ4
                        T4 = (Q4-Q3) / (P4-P3)
                        P14 = P4
                        Q14 = Q4
                    end if
                    Gyp = P14
                case(7) !Stage 2' (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P5 = LD
                    Q5 = LS
                    T1 = IT1
                    P6 = P5 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    Q6 = Q5 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    T6 = T1 !�����̃X�e�[�W2�͈ꎞ�I�ɏ�������B
                    P7 = P7
                    Q7 = Q7
                    T7 = T7
                    P8 = P8
                    Q8 = Q8
                    T8 = T8
                    P18 = P18
                    Q18 = Q18
                    Gyn = Gyn
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P5 + (Q1-Q5)/T1
                    P2 = IP2
                    Q2 = IQ2
                    T2 = (Q2-Q1)/(P2-P1)
                    Q3 = IQ3
                    P3 = IP3
                    T3 = IT3
                    if     (Pmax <  IP4) then !�����̋Ȃ��~����o��, Rule 1
                        P4 = IP4
                        Q4 = IQ4
                        T4 = IT4
                        P14 = P4
                        Q14 = Q4
                    else if(Pmax >= IP4) then !�����̋Ȃ��~���@�o��, Rule 3
                        T4 = FuncAlpha(IP2,IP4,Pmax) * IT4
                        Q4 = IQ4
                        P4 = P3 + (Q4-Q3)/T4
                        P14 = Pmax                    !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                        Q14 = Q4 + IT2 * (P14-P4) !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                    end if
                    Gyp = P14
                case(8) !Stage 3 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P5 = LD
                    Q5 = LS
                    T1 = IT1
                    P6 = P5 !������Stage2�͈ꎞ�I�ɏ�������
                    Q6 = Q5 !������Stage2�͈ꎞ�I�ɏ�������
                    T6 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P7 = P5 !������Stage2'�͈ꎞ�I�ɏ�������
                    Q7 = Q5 !������Stage2'�͈ꎞ�I�ɏ�������
                    T7 = T1 !������Stage2'�͈ꎞ�I�ɏ�������
                    P8 = P8
                    Q8 = Q8
                    T8 = T8
                    P18 = P18
                    Q18 = Q18
                    Gyn = Gyn
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P5 + (Q1-Q5)/T1
                    P2 = IP2
                    Q2 = IQ2
                    T2 = (Q2-Q1)/(P2-P1)
                    Q3 = IQ3
                    P3 = IP3
                    T3 = IT3
                    if     (Pmax <  IP4) then !�����̋Ȃ��~����o��, Rule 1
                        P4 = IP4
                        Q4 = IQ4
                        T4 = IT4
                        P14 = P4
                        Q14 = Q4
                    else if(Pmax >= IP4) then !�����̋Ȃ��~���@�o��, Rule 3
                        T4 = FuncAlpha(IP2,IP4,Pmax) * IT4
                        Q4 = IQ4
                        P4 = P3 + (Q4-Q3)/T4
                        P14 = Pmax                    !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                        Q14 = Q4 + IT2 * (P14-P4) !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                    end if
                    Gyp = P14
                case(9) !Stage 4 (�׏d�������̃T�C�h)
                    !----�����̍��i�Ȑ��͍X�V�Ȃ�----!
                    P5 = LD
                    Q5 = LS
                    T1 = IT1
                    P6 = P5 !������Stage2�͈ꎞ�I�ɏ�������
                    Q6 = Q6 !������Stage2�͈ꎞ�I�ɏ�������
                    T6 = T1 !������Stage2�͈ꎞ�I�ɏ�������
                    P7 = P5 !������Stage2'�͈ꎞ�I�ɏ�������
                    Q7 = Q5 !������Stage2'�͈ꎞ�I�ɏ�������
                    T7 = T1 !������Stage2'�͈ꎞ�I�ɏ�������
                    P8 = P5 !������Stage3�͈ꎞ�I�ɏ�������
                    Q8 = Q5 !������Stage3�͈ꎞ�I�ɏ�������
                    T8 = T1 !������Stage3�͈ꎞ�I�ɏ�������
                    P18 = P8
                    Q18 = Q8
                    Gyn = IP8 !������Stage4�̕ό`�J�n�_�͏����Ȃ��~���_�Ɠ���
                    !----�����̍��i�Ȑ��͍X�V����----!
                    Q1 = IQ1
                    P1 = P5 + (Q1-Q5)/T1
                    P2 = IP2
                    Q2 = IQ2
                    T2 = (Q2-Q1)/(P2-P1)
                    Q3 = IQ3
                    P3 = IP3
                    T3 = IT3
                    if     (Pmax <  IP4) then !�����̋Ȃ��~����o��, Rule 1
                        P4 = IP4
                        Q4 = IQ4
                        T4 = IT4
                        P14 = P4
                        Q14 = Q4
                    else if(Pmax >= IP4) then !�����̋Ȃ��~���@�o��, Rule 3
                        T4 = FuncAlpha(IP2,IP4,Pmax) * IT4
                        Q4 = IQ4
                        P4 = P3 + (Q4-Q3)/T4
                        P14 = Pmax                    !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                        Q14 = Q4 + IT2 * (P14-P4) !�����̋Ȃ��~���_��4����14�ֈړ�����(4����14�̊Ԃ����镔��)
                    end if
                    Gyp = P14
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
    if     ( D >= P14            ) then ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
        stage = 5
    else if( D >= P3            ) then  ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
        stage = 4
    else if( D >= P2            ) then  ! ������Stage2'(�A���J�[�{���g��BPL�ƐڐG���ĕ~�v���[�g�ƃ����^��������o���܂ł̖��C��)
        stage = 3
    else if( D >= P1            ) then  ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
        stage = 2
    else if( D >  P5.and.D < P1 ) then  ! Stage1 (����o���܂�)
        stage = 1
    else if( D <= P18            ) then ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
        stage = 9
    else if( D <= P7            ) then  ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
        stage = 8
    else if( D <= P6            ) then  ! ������Stage2'(�A���J�[�{���g��BPL�ƐڐG���ĕ~�v���[�g�ƃ����^��������o���܂ł̖��C��)
        stage = 7
    else if( D <= P5            ) then  ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
        stage = 6
    else
        print *,"Stage error 9999"
        print *,"D :" ,D
        print *,"S :" ,S
        print *,"T :" ,T
        print *,"P1:",P1
        print *,"P2:",P2
        print *,"P3:",P3
        print *,"P4:",P4
        print *,"P14:",P14
        print *,"P5:",P5
        print *,"P6:",P6
        print *,"P7:",P7
        print *,"P8:",P8
        print *,"P18:",P18
        print *,"Q1:",Q1
        print *,"Q2:",Q2
        print *,"Q3:",Q3
        print *,"Q4:",Q4
        print *,"Q14:",Q14
        print *,"Q5:",Q5
        print *,"Q6:",Q6
        print *,"Q7:",Q7
        print *,"Q8:",Q8
        print *,"Q18:",Q18
        print *,"T1:",T1
        print *,"T2:",T2
        print *,"T3:",T3
        print *,"T4:",T4
        print *,"T6:",T6
        print *,"T7:",T7
        print *,"T8:",T8
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
        case(5) ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
            S = AnchorGeoNonlinear(nn,mu,IP3,Gyp,D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = IT2
        case(4) ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
            if( D <  P14.and.D >= P4) then  ! ������Stage4+Rule3 (�A���J�[�{���g���ēx��������܂ł̃X���b�v��)
                call LinearSlope(S,T,D,IT2,P4,Q4)
            else
                call LinearSlope(S,T,D,T4,P3,Q3)
            end if
        case(3) ! ������Stage2'(�A���J�[�{���g��BPL�ƐڐG���ĕ~�v���[�g�ƃ����^��������o���܂ł̖��C��)
            call LinearSlope(S,T,D,T3,P2,Q2)
        case(2) ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
            call LinearSlope(S,T,D,T2,P1,Q1)
        case(1) ! Stage1 (����o���܂�)
            call LinearSlope(S,T,D,T1,P1,Q1)
        case(6) ! ������Stage2�iBPL������o���ăA���J�[�{���g�ƐڐG����܂�)
            call LinearSlope(S,T,D,T6,P5,Q5)
        case(7) ! ������Stage2'(�A���J�[�{���g��BPL�ƐڐG���ĕ~�v���[�g�ƃ����^��������o���܂ł̖��C��)
            call LinearSlope(S,T,D,T7,P6,Q6)
        case(8) ! ������Stage3 (�A���J�[�{���g��BPL�ƐڐG���ċȂ��~������܂ł̑ϗ͏㏸��)
            if     ( D <= P8.and.D > P18 ) then  ! ������Stage4+Rule3 (�A���J�[�{���g���ēx��������܂ł̃X���b�v��)
                call LinearSlope(S,T,D,IT2,P8,Q8)
            else
                call LinearSlope(S,T,D,T8,P7,Q7)
            end if
        case(9) ! ������Stage4 (�A���J�[�{���g�̈����ɔ����􉽔���`���ɂ��ϗ͏㏸��)
            S = -1.0d0 * AnchorGeoNonlinear(nn,mu,IP3,-1.0d0*Gyn,-1.0d0*D,C4,C1,C2,C3,C5)
            if(dabs(D-LD) /= 0.0d0) T = dabs(S-LS)/dabs(D-LD)
            if(T <= 0.0d0) T = IT2
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
    if(Ms < 0 .or. Ns >= CC5) then !�Y���֌W���ɂ�鐧��
        !print *, Ms,Ns
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