#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	vec3 root_guideV = m_pSkeleton->getRootNode()->getGlobalTranslation();
	vec3 root_worldV = m_Guide.getLocal2Global() * root_guideV;

	vec3 guide_worldV = m_Guide.getGlobalTranslation();
	vec3 tar = (guideTargetPos - guide_worldV).Normalize();

	m_Guide.setGlobalTranslation(vec3(root_worldV[0], 0.f, root_worldV[2]));

	m_Guide.setGlobalRotation(mat3(vec3(0.f, 1.f, 0.f).Cross(tar), vec3(0.f,1.f,0.f), tar).Transpose());
	m_pSkeleton->update();
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint* root = m_pSkeleton->getRootNode();
	vec3 hL = root->getLocalTranslation() + vec3(0, leftHeight, 0);
	vec3 hR = root->getLocalTranslation() + vec3(0, rightHeight, 0);
	root->setLocalTranslation(Max(hL, hR));


	// 2.	Update the character with Limb-based IK 
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal

		ATarget lTarget = ATarget();
		vec3 xAxis = (leftFoot->getLocalRotation().GetCol(0)).Normalize();
		vec3 newZ = (leftNormal.Cross(xAxis)).Normalize();
		lTarget.setLocalRotation(mat3(xAxis, leftNormal.Normalize(), newZ).Transpose());
		
		
		m_IKController->computeLimbIK(lTarget, m_IKController->createIKchain(leftFoot->getID(), 3, m_pSkeleton), 
			xAxis, m_pSkeleton);
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATarget rTarget = ATarget();
		vec3 xAxis = (rightFoot->getLocalRotation().GetCol(0)).Normalize();
		vec3 newZ = (rightNormal.Cross(xAxis)).Normalize();
		rTarget.setLocalRotation(mat3(xAxis, rightNormal.Normalize(), newZ).Transpose());
		
		m_IKController->computeLimbIK(rTarget, m_IKController->createIKchain(rightFoot->getID(), 3, m_pSkeleton),
			xAxis, m_pSkeleton);
	}
	m_pSkeleton->update();
}
