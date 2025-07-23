-- Author: Kyle Damm
-- Date: 21-September-2021
--
-- Transcribed from page 4 of
-- An analysis of combustion studies in shock expansion tunnels and reflected shock tunnels
-- Jachimowski, 1992

Config{
   odeStep = {method='alpha-qss'},
   tightTempCoupling = true,
}

Reaction{
   'H + OH + M <=> H2O + M',
   fr={'Arrhenius', A=8.62e21, n=-2.0, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r6'
}

Reaction{
   'H + H + M <=> H2 + M',
   fr={'Arrhenius', A=7.3e17, n=-1.0, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r7'
}

Reaction{
   'H + O + M <=> OH + M',
   fr={'Arrhenius', A=2.6e16, n=-0.6, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r8'
}

Reaction{
   'O + O + M <=> O2 + M',
   fr={'Arrhenius', A=1.1e17, n=-1.0, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r9'
}

Reaction{
   'H + O2 + M <=> HO2 + M',
   fr={'Arrhenius', A=2.3e18, n=-1.0, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r10'
}
Reaction{
   'OH + OH + M <=> H2O2 + M',
   fr={'Arrhenius', A=1.6e22, n=-2.0, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r19'
}

Reaction{
   'N + N + M <=> N2 + M',
   fr={'Arrhenius', A=2.8e17, n=-0.8, C=0.0},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r20'
}

Reaction{
   'H + NO + M <=> HNO + M',
   fr={'Arrhenius', A=5.4e15, n=0.0, C=-3.02133e02},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r24'
}

Reaction{
   'M + NO2 <=> NO + O + M',
   fr={'Arrhenius', A=1.16e16, n=0.0, C=3.32347e04},
   efficiencies={['H2']=2.5,['H2O']=16.0},
   label='r33'
}

Reaction{
   'H2 + O2 <=> H + HO2',
   fr={'Arrhenius', A=7.0e13, n=0.0, C=2.86020e04},
   label='r1'
}

Reaction{
   'H + O2 <=> OH + O',
   fr={'Arrhenius', A=2.2e14, n=0.0, C=8.45974e03},
   label='r2'
}

Reaction{
   'O + H2 <=> OH + H',
   fr={'Arrhenius', A=5.06e4, n=2.67, C=3.16737e03},
   label='r3'
}

Reaction{
   'OH + H2 <=> H2O + H',
   fr={'Arrhenius', A=2.16e08, n=1.51, C=1.72720e03},
   label='r4'
}

Reaction{
   'OH + OH <=> H2O + O',
   fr={'Arrhenius', A=1.5e09, n=1.14, C=0.0},
   label='r5'
}


Reaction{
   'HO2 + H <=> OH + OH',
   fr={'Arrhenius', A=1.5e14, n=0.0, C=5.03556e02},
   label='r11'
}

Reaction{
   'HO2 + O <=> O2 + OH',
   fr={'Arrhenius', A=2.0e13, n=0.0, C=0.0},
   label='r12'
}

Reaction{
   'HO2 + OH <=> H2O + O2',
   fr={'Arrhenius', A=2.0e13, n=0.0, C=0.0},
   label='r13'
}

Reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=2.0e12, n=0.0, C=0.0},
   label='r14'
}

Reaction{
   'H + H2O2 <=> H2 + HO2',
   fr={'Arrhenius', A=1.7e12, n=0.0, C=1.90344e03},
   label='r15'
}

Reaction{
   'H + H2O2 <=> OH + H2O',
   fr={'Arrhenius', A=1.0e13, n=0.0, C=1.80273e03},
   label='r16'
}

Reaction{
   'O + H2O2 <=> OH + HO2',
   fr={'Arrhenius', A=2.8e13, n=0.0, C=3.22276e03},
   label='r17'
}

Reaction{
   'OH + H2O2 <=> H2O + HO2',
   fr={'Arrhenius', A=7.0e12, n=0.0, C=7.22602e02},
   label='r18'
}


Reaction{
   'N + O2 <=> NO + O',
   fr={'Arrhenius', A=6.4e09, n=1.0, C=3.17240e03},
   label='r21'
}

Reaction{
   'N + NO <=> N2 + O',
   fr={'Arrhenius', A=1.6e13, n=0.0, C=0.0},
   label='r22'
}

Reaction{
   'N + OH <=> NO + H',
   fr={'Arrhenius', A=6.3e11, n=0.5, C=0.0},
   label='r23'
}

Reaction{
   'H + HNO <=> NO + H2',
   fr={'Arrhenius', A=4.8e12, n=0.0, C=0.0},
   label='r25'
}

Reaction{
   'O + HNO <=> NO + OH',
   fr={'Arrhenius', A=5.0e11, n=0.5, C=0.0},
   label='r26'
}

Reaction{
   'OH + HNO <=> NO + H2O',
   fr={'Arrhenius', A=3.6e13, n=0.0, C=0.0},
   label='r27'
}

Reaction{
   'HO2 + HNO <=> NO + H2O2',
   fr={'Arrhenius', A=2.0e12, n=0.0, C=0.0},
   label='r28'
}

Reaction{
   'HO2 + NO <=> NO2 + OH',
   fr={'Arrhenius', A=3.4e12, n=0.0, C=-1.30924e02},
   label='r29'
}

Reaction{
   'HO2 + NO <=> HNO + O2',
   fr={'Arrhenius', A=2.0e11, n=0.0, C=5.03556e02},
   label='r30'
}

Reaction{
   'H + NO2 <=> NO + OH',
   fr={'Arrhenius', A=3.5e14, n=0.0, C=7.55334e02},
   label='r31'
}

Reaction{
   'O + NO2 <=> NO + O2',
   fr={'Arrhenius', A=1.0e13, n=0.0, C=3.02133e02},
   label='r32'
}
