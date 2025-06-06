#pragma once
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/VectorUtil.h>

namespace WWTools
{
    using namespace ROOT::Math;

    // takes W and fermion 4 vector in the WW/CMS frame and input coordinate system specified by y and z
    template <class CoordSystem>
    LorentzVector<CoordSystem>
    starVectorBasis(const LorentzVector<CoordSystem>& W_lvec, const LorentzVector<CoordSystem>& lep_lvec, const XYZVector& y, const XYZVector& z)
    {
        // the obvious third dimension
        const auto x = y.Cross(z);

        // .Inverse() to obtain passive rotation (rotate coord system not vector)
        const auto rot = Rotation3D(x, y, z).Inverse();

        // get the boost to boost something into the W CMS, boost and then rotate
        const auto W_boost = W_lvec.BoostToCM();
        const auto boosted_lvec = VectorUtil::boost(lep_lvec, W_boost);
        const auto transformed_lvec = rot(boosted_lvec);

        return transformed_lvec;
    }

    template <class CoordSystem>
    LorentzVector<CoordSystem>
    starVectorHagiwara(const LorentzVector<CoordSystem>& W_lvec, const LorentzVector<CoordSystem>& lep_lvec, const XYZVector& e_minus_vec, const LorentzVector<CoordSystem>& W_minus_lvec)
    {
        // z* points along the W flight direction in the lab frame
        const auto z = W_minus_lvec.Vect().Unit();
        // y* points in the direction of e- x W
        const auto y = e_minus_vec.Cross(z).Unit();

        return starVectorBasis(W_lvec, lep_lvec, y, z);
    }

    template <class CoordSystem>
    LorentzVector<CoordSystem>
    starVector(const LorentzVector<CoordSystem>& W_lvec, const LorentzVector<CoordSystem>& lep_lvec, const XYZVector& e_minus_vec = XYZVector(0, 0, 1))
    {
        // z* points along the W flight direction in the lab frame
        const auto z = W_lvec.Vect().Unit();
        // y* points in the direction of e- x W
        const auto y = e_minus_vec.Cross(z).Unit();

        return starVectorBasis(W_lvec, lep_lvec, y, z);
    }

    // XXX: don't use, this is certainly wrong!
    using namespace ROOT::Math;
    template <class CoordSystem>
    LorentzVector<CoordSystem>
    starVectorOO(const LorentzVector<CoordSystem>& W_lvec, const LorentzVector<CoordSystem>& lep_lvec, const XYZVector& e_minus_vec = XYZVector(0, 0, 1))
    {
        // z* points along the W flight direction in the lab frame
        // const auto z = W_lvec.Vect().Unit();
        // y* points in the direction of e- x W
        // const auto y = e_minus_vec.Cross(z).Unit();
        // the obvious third dimension
        // const auto x = y.Cross(z);

        // don't rotate
        // .Inverse() to obtain passive rotation (rotate coord system not vector)
        // const auto rot = Rotation3D(x, y, z).Inverse();

        // get the boost to boost something into the W CMS, boost and then rotate
        const auto W_boost = W_lvec.BoostToCM();
        const auto boosted_lvec = VectorUtil::boost(lep_lvec, W_boost);
        // const auto transformed_lvec = rot(boosted_lvec);

        // return transformed_lvec;
        return boosted_lvec;
    }
}