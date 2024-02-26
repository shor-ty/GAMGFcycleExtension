/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "GAMGSolverNew.H"
#include "SubField.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::GAMGSolverNew::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // Calculate A.psi used to calculate the initial residual
    solveScalarField Apsi(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    solveScalarField finestCorrection(psi.size());

    // Calculate normalisation factor
    solveScalar normFactor =
        this->normFactor(psi, tsource(), Apsi, finestCorrection);

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    solveScalarField finestResidual(tsource() - Apsi);

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        true
    );

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag
    (
        finestResidual,
        matrix().mesh().comm()
    )/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // Create coarse grid correction fields
        PtrList<solveScalarField> coarseCorrFields;

        // Create coarse grid sources
        PtrList<solveScalarField> coarseSources;

        // Create the smoothers for all levels
        PtrList<lduMatrix::smoother> smoothers;

        // Scratch fields if processor-agglomerated coarse level meshes
        // are bigger than original. Usually not needed
        solveScalarField scratch1;
        solveScalarField scratch2;

        // Initialise the above data structures
        initVcycle
        (
            coarseCorrFields,
            coarseSources,
            smoothers,
            scratch1,
            scratch2
        );

        do
        {
            //Info<< "Calculate GAMG Cycle" << nl;
            if (cycleMode_ == "VCycle")
            {
                Vcycle
                (
                    smoothers,
                    psi,
                    source,
                    Apsi,
                    finestCorrection,
                    finestResidual,

                    (scratch1.size() ? scratch1 : Apsi),
                    (scratch2.size() ? scratch2 : finestCorrection),

                    coarseCorrFields,
                    coarseSources,
                    cmpt
                );
            } 
            else if (cycleMode_ == "FCycle")
            {
                FCycle
                (
                    smoothers,
                    psi,
                    source,
                    Apsi,
                    finestCorrection,
                    finestResidual,

                    (scratch1.size() ? scratch1 : Apsi),
                    (scratch2.size() ? scratch2 : finestCorrection),

                    coarseCorrFields,
                    coarseSources,
                    cmpt
                );
            }
            else
            {
                std::terminate();
            }

            // Calculate finest level residual field
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            finestResidual = tsource();
            finestResidual -= Apsi;

            solverPerf.finalResidual() = gSumMag
            (
                finestResidual,
                matrix().mesh().comm()
            )/normFactor;

            if ((log_ >= 2) || (debug >= 2))
            {
                solverPerf.print(Info.masterStream(matrix().mesh().comm()));
            }
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        false
    );

    return solverPerf;
}


void Foam::GAMGSolverNew::restriction
(
    const label leveli,
    const label coarsestLevel,
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Restrict from level " << leveli + 1 << " to " << leveli + 2 << endl;
    }

    // If the optional pre-smoothing sweeps are selected
    // smooth the coarse-grid field for the restricted source
    if (nPreSweeps_)
    {
        coarseCorrFields[leveli] = 0.0;

        smoothers[leveli + 1].scalarSmooth
        (
            coarseCorrFields[leveli],
            coarseSources[leveli],  //coarseSource,
            cmpt,
            min
            (
                nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                maxPreSweeps_
            )
        );

        solveScalarField::subField ACf
        (
            scratch1,
            coarseCorrFields[leveli].size()
        );

        // Scale coarse-grid correction field
        // but not on the coarsest level because it evaluates to 1
        if (scaleCorrection_ && leveli < coarsestLevel - 1)
        {
            scale
            (
                coarseCorrFields[leveli],
                const_cast<solveScalarField&>
                (
                    ACf.operator const solveScalarField&()
                ),
                matrixLevels_[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                coarseSources[leveli],
                cmpt
            );
        }

        // Correct the residual with the new solution
        matrixLevels_[leveli].Amul
        (
            const_cast<solveScalarField&>
            (
                ACf.operator const solveScalarField&()
            ),
            coarseCorrFields[leveli],
            interfaceLevelsBouCoeffs_[leveli],
            interfaceLevels_[leveli],
            cmpt
        );

        coarseSources[leveli] -= ACf;
    }

    // Residual is equal to source
    agglomeration_.restrictField
    (
        coarseSources[leveli + 1],
        coarseSources[leveli],
        leveli + 1,
        true
    );
}


void Foam::GAMGSolverNew::prolongation
(
    const label leveli,
    const label coarsestLevel,
    solveScalarField& dummyField,
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    solveScalarField& scratch2,
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Prolongate from level " << leveli + 2 << " to " << leveli + 1 << " Scaling = ";
    }

    // Create a field for the pre-smoothed correction field
    // as a sub-field of the finestCorrection which is not
    // currently being used
    solveScalarField::subField preSmoothedCoarseCorrField
    (
        scratch2,
        coarseCorrFields[leveli].size()
    );

    // Only store the preSmoothedCoarseCorrField if pre-smoothing is
    // used
    if (nPreSweeps_)
    {
        preSmoothedCoarseCorrField = coarseCorrFields[leveli];
    }

    agglomeration_.prolongField
    (
        coarseCorrFields[leveli],
        (
            coarseCorrFields.set(leveli + 1)
          ? coarseCorrFields[leveli + 1]
          : dummyField              // dummy value
        ),
        leveli + 1,
        true
    );


    // Create A.psi for this coarse level as a sub-field of Apsi
    solveScalarField::subField ACf
    (
        scratch1,
        coarseCorrFields[leveli].size()
    );
 
    solveScalarField& ACfRef =
        const_cast
        <
            solveScalarField&
        >(ACf.operator const solveScalarField&());

    if (interpolateCorrection_) //&& leveli < coarsestLevel - 2)
    {
        if (coarseCorrFields.set(leveli+1))
        {
            interpolate
            (
                coarseCorrFields[leveli],
                ACfRef,
                matrixLevels_[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                agglomeration_.restrictAddressing(leveli + 1),
                coarseCorrFields[leveli + 1],
                cmpt
            );
        }
        else
        {
            interpolate
            (
                coarseCorrFields[leveli],
                ACfRef,
                matrixLevels_[leveli],
                interfaceLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                cmpt
            );
        }
    }

    // Scale coarse-grid correction field
    // but not on the coarsest level because it evaluates to 1
    if
    (
        scaleCorrection_
     && (interpolateCorrection_ || leveli < coarsestLevel - 1)
    )
    {
        scale
        (
            coarseCorrFields[leveli],
            ACfRef,
            matrixLevels_[leveli],
            interfaceLevelsBouCoeffs_[leveli],
            interfaceLevels_[leveli],
            coarseSources[leveli],
            cmpt
        );
    }

    // Only add the preSmoothedCoarseCorrField if pre-smoothing is
    // used
    if (nPreSweeps_)
    {
        coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
    }

    smoothers[leveli + 1].scalarSmooth
    (
        coarseCorrFields[leveli],
        coarseSources[leveli],  //coarseSource,
        cmpt,
        min
        (
            nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
            maxPostSweeps_
        )
    );

    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< endl;
    }
}


void Foam::GAMGSolverNew::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& psi,
    const scalarField& source,
    solveScalarField& Apsi,
    solveScalarField& finestCorrection,
    solveScalarField& finestResidual,

    solveScalarField& scratch1,
    solveScalarField& scratch2,

    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Restrict from level 0 to 1" << endl;
    }

    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli <= coarsestLevel; leveli++)
    {
        if (coarseSources.set(leveli + 1))
        {
            restriction
            (
                leveli,
                coarsestLevel,
                smoothers,
                scratch1,
                coarseCorrFields,
                coarseSources,
                cmpt
            );
        }
    }

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    //Info<< "+" << endl;
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    solveScalarField dummyField(0);

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        if (coarseCorrFields.set(leveli))
        {
            prolongation
            (
                leveli,
                coarsestLevel,
                dummyField,
                smoothers,
                scratch1,
                scratch2,
                coarseCorrFields,
                coarseSources,
                cmpt
            );
        }
    }

    // Prolong the finest level correction
    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Prolongate from level 1 to 0" << endl;
    }

    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0],
            cmpt
        );
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
}


void Foam::GAMGSolverNew::FCycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& psi,
    const scalarField& source,
    solveScalarField& Apsi,
    solveScalarField& finestCorrection,
    solveScalarField& finestResidual,

    solveScalarField& scratch1,
    solveScalarField& scratch2,

    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;
    if  ((log_ >= 2) || (debug >= 2))
    {
        Info<< "Restrict from level 0 to 1" << endl;
    }

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    // --------------------------------------------------------------
    // Left side
    // --------------------------------------------------------------
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        if (coarseSources.set(leveli + 1))
        {
            restriction
            (
                leveli,
                coarsestLevel,
                smoothers,
                scratch1,
                coarseCorrFields,
                coarseSources,
                cmpt
            );
        }
    }

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< endl;
    }

    // Solve Coarsest level with either an iterative or direct solver
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // --------------------------------------------------------------
    // Inbetween
    // --------------------------------------------------------------
    // Go up one level and back to coarsestLevel
    // Go up two levels and back to coarsestLevel
    // Go up three levels and back to coarsestLevel
    // ...
    // Until we are on the right side
    // --------------------------------------------------------------
    for 
    (
        label toLeveli = coarsestLevel - 1;
        toLeveli >= 1;
        toLeveli--
    ) 
    {
        solveScalarField dummyField(0);
        // Prolongate to level toLeveli
        for 
        (
            label prolongLeveli = coarsestLevel - 1;
            prolongLeveli >= toLeveli;
            prolongLeveli--
        ) 
        {
            if (coarseCorrFields.set(prolongLeveli))
            {
                prolongation
                (
                    prolongLeveli,
                    coarsestLevel,
                    dummyField,
                    smoothers,
                    scratch1,
                    scratch2,
                    coarseCorrFields,
                    coarseSources,
                    cmpt
                );
            }
        }

        {

            if ((log_ >= 2) || (debug >= 2))
            {
                Pout<< "Update residuals on level " << toLeveli << endl;
            }

            solveScalarField::subField ACf
            (
                scratch1,
                coarseCorrFields[toLeveli].size()
            );

            // Correct the residual with the new solution
            matrixLevels_[toLeveli].Amul
            (
                const_cast<solveScalarField&>
                (
                    ACf.operator const solveScalarField&()
                ),
                coarseCorrFields[toLeveli],
                interfaceLevelsBouCoeffs_[toLeveli],
                interfaceLevels_[toLeveli],
                cmpt
            );

            coarseSources[toLeveli] -= ACf;
        }

        // Restrict to coarsestLevel again and solve on coarsestLevel
        for 
        (
            label restrictLeveli = toLeveli;
            restrictLeveli <= coarsestLevel - 1;
            restrictLeveli++
        ) 
        {
            if (coarseCorrFields.set(restrictLeveli))
            {
                restriction 
                (
                    restrictLeveli,
                    coarsestLevel,
                    smoothers,
                    scratch1,
                    coarseCorrFields,
                    coarseSources,
                    cmpt
                );
            }
        }

        // Solve Coarsest level with either an iterative or direct solver
        if (coarseCorrFields.set(coarsestLevel))
        {
            solveCoarsestLevel
            (
                coarseCorrFields[coarsestLevel],
                coarseSources[coarsestLevel]
            );
        }
    }

    // --------------------------------------------------------------
    // Right side
    // --------------------------------------------------------------
    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)
    solveScalarField dummyField(0);

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        if (coarseCorrFields.set(leveli))
        {
            prolongation
            (
                leveli,
                coarsestLevel,
                dummyField,
                smoothers,
                scratch1,
                scratch2,
                coarseCorrFields,
                coarseSources,
                cmpt
            );
        }
    }

    // Prolong the finest level correction
    if  ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Prolongate from level 1 to 0" << endl;
    }

    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0],
            cmpt
        );
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
}


void Foam::GAMGSolverNew::initVcycle
(
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    solveScalarField& scratch2
) const
{
    label maxSize = matrix_.diag().size();

    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Create the smoother for the finest level
    smoothers.set
    (
        0,
        lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        )
    );

    forAll(matrixLevels_, leveli)
    {
        if (agglomeration_.nCells(leveli) >= 0)
        {
            label nCoarseCells = agglomeration_.nCells(leveli);

            coarseSources.set(leveli, new solveScalarField(nCoarseCells));
        }

        if (matrixLevels_.set(leveli))
        {
            const lduMatrix& mat = matrixLevels_[leveli];

            label nCoarseCells = mat.diag().size();

            maxSize = max(maxSize, nCoarseCells);

            coarseCorrFields.set(leveli, new solveScalarField(nCoarseCells));

            smoothers.set
            (
                leveli + 1,
                lduMatrix::smoother::New
                (
                    fieldName_,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevelsIntCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    controlDict_
                )
            );
        }
    }

    if (maxSize > matrix_.diag().size())
    {
        // Allocate some scratch storage
        scratch1.setSize(maxSize);
        scratch2.setSize(maxSize);
    }
}


Foam::dictionary Foam::GAMGSolverNew::PCGsolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PCG; preconditioner DIC;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::GAMGSolverNew::PBiCGStabSolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PBiCGStab; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


void Foam::GAMGSolverNew::solveCoarsestLevel
(
    solveScalarField& coarsestCorrField,
    const solveScalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    const label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();

    if (directSolveCoarsest_)
    {
        PrecisionAdaptor<scalar, solveScalar> tcorrField(coarsestCorrField);

        coarsestLUMatrixPtr_->solve
        (
            tcorrField.ref(),
            ConstPrecisionAdaptor<scalar, solveScalar>(coarsestSource)()
        );
    }
    //else if
    //(
    //    agglomeration_.processorAgglomerate()
    // && procMatrixLevels_.set(coarsestLevel)
    //)
    //{
    //    //const labelList& agglomProcIDs = agglomeration_.agglomProcIDs
    //    //(
    //    //    coarsestLevel
    //    //);
    //    //
    //    //scalarField allSource;
    //    //
    //    //globalIndex cellOffsets;
    //    //if (Pstream::myProcNo(coarseComm) == agglomProcIDs[0])
    //    //{
    //    //    cellOffsets.offsets() =
    //    //        agglomeration_.cellOffsets(coarsestLevel);
    //    //}
    //    //
    //    //cellOffsets.gather
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    coarsestSource,
    //    //    allSource
    //    //);
    //    //
    //    //scalarField allCorrField;
    //    //solverPerformance coarseSolverPerf;
    //
    //    label solveComm = agglomeration_.procCommunicator(coarsestLevel);
    //
    //    coarsestCorrField = 0;
    //    solverPerformance coarseSolverPerf;
    //
    //    if (Pstream::myProcNo(solveComm) != -1)
    //    {
    //        const lduMatrix& allMatrix = procMatrixLevels_[coarsestLevel];
    //
    //        {
    //            Pout<< "** Master:Solving on comm:" << solveComm
    //                << " with procs:" << UPstream::procID(solveComm) << endl;
    //
    //            if (allMatrix.asymmetric())
    //            {
    //                coarseSolverPerf = PBiCGStab
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PBiCGStabSolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //            else
    //            {
    //                coarseSolverPerf = PCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PCGsolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //        }
    //    }
    //
    //    Pout<< "done master solve." << endl;
    //
    //    //// Scatter to all processors
    //    //coarsestCorrField.setSize(coarsestSource.size());
    //    //cellOffsets.scatter
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    allCorrField,
    //    //    coarsestCorrField
    //    //);
    //
    //    if (debug >= 2)
    //    {
    //        coarseSolverPerf.print(Info.masterStream(coarseComm));
    //    }
    //
    //    Pout<< "procAgglom: coarsestSource   :" << coarsestSource << endl;
    //    Pout<< "procAgglom: coarsestCorrField:" << coarsestCorrField << endl;
    //}
    else
    {
        coarsestCorrField = 0;
        const solverPerformance coarseSolverPerf
        (
            coarsestSolverPtr_->scalarSolve
            (
                coarsestCorrField,
                coarsestSource
            )
        );

        if ((log_ >= 2) || debug)
        {
            coarseSolverPerf.print(Info.masterStream(coarseComm));
        }
    }
}


// ************************************************************************* //
