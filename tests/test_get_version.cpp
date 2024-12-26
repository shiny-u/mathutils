//
// Created by frongere on 06/04/2021.
//

#include <iostream>
#include <MathUtils/version.h>

#include <gtest/gtest.h>


using namespace mathutils::git;


TEST(GetVersion, Demo) {
    std::cout << ProjectName() << " Version: " << GetNormalizedVersionString() << std::endl;
    std::cout << "Has Uncommitted changes: " << AnyUncommittedChanges() << std::endl;
    std::cout << "Commit: " << CommitSHA1() << std::endl;
    std::cout << "Date: " << CommitDate() << std::endl;
    std::cout << "Author: " << AuthorName() << std::endl;
}
